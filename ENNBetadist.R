################################################################################
####                                                                        ####
####          New deep learning-based methods                               ####
####        for visualizing ecosystem properties                            ####
####       using environmental DNA metabarcoding data                       ####
####                                                                        ####
####                                                                        ####
####                                ENNBetaDist                             ####
####                                                                        ####
####   Letizia Lamperti et al. 2023                                         ####
####                                                                        ####
####   letizia.lamperti@ephe.psl.eu                                         ####
################################################################################


### Libraries loading
library(keras)
library(tensorflow)
library(MLmetrics)
library(gdata)
library(dplyr)

### Set working directory
setwd("")

#define conda environment
use_condaenv("")

#importing KERAS in R
library(keras)
K <- keras::backend()

#import data from Telonemia or Dictyochophyceae fold 
Database_W5_OTU_occurences_1 <- read.csv("Data/Telonemia/Database_W5_OTU_occurences_Telonemia.csv")
Data <- as.data.frame(Database_W5_OTU_occurences_1[,4:337])

Data <- as.data.frame(t(Data))
Data <- subset(Data, rowSums(Data)>0) #eliminating empty-samples 
n_MOTU <- dim(Data)[2]
n_samples <- dim(Data)[1]

#normalizing data 
Data <- Data / max(Data)
dim(Data)

#import Genetic info 
genetic_info <- read.csv("Data/Telonemia/genetic_info_1_bs_100_Telonemia.csv")
genetic_info <- as.matrix(genetic_info)
#normalizing data 
genetic_info <- genetic_info - min(genetic_info)
genetic_info <- genetic_info / max(genetic_info)
genetic_info <- as.data.frame(genetic_info)

#compute beta distances between points 
library(betapart)
Data_presence <- as.matrix((Data>0)+0)
jacc_dist <- beta.pair(Data_presence, index.family="jaccard")
matrix_jacc_dist <- as.data.frame(as.numeric(jacc_dist$beta.jac))

beta_values <- as.matrix(matrix_jacc_dist)


#Writing input of the NN
#input 1: i-th row repeated (dim(samples) - i) times
#input 2: Data - i rows

#INPUT 1
n_samples <- dim(Data)[1]
input_1 <- NULL

library(vctrs)
for (i in 1: n_samples) {
  sample <- as.vector(Data[i,])
  repeated_sample <-vec_rep(sample , n_samples - i)
  input_1 <- rbind(input_1, repeated_sample)
}

m_input_1 <- as.matrix(input_1[, 1:n_MOTU])
dim(m_input_1)

#INPUT 2
input_2 <- NULL

for (i in 1: n_samples) {
  Data_sample <- NULL
  Data_sample <- as.vector(Data[-(1:i),])
  input_2 <- rbind(input_2, Data_sample)
}

m_input_2 <- as.matrix(input_2[, 1:n_MOTU])
dim(m_input_2)

#GENETIC INPUT 1
n_samples <- dim(genetic_info)[1] 
genetic_info_1 <- NULL

library(vctrs)
for (i in 1: n_samples) {
  sample <- as.vector(genetic_info[i,])
  repeated_sample <-vec_rep(sample , n_samples - i)
  genetic_info_1 <- rbind(genetic_info_1, repeated_sample)
}

genetic_info_1 <- as.matrix(genetic_info_1)
dim(genetic_info_1)

#GENETIC INPUT 2
genetic_info_2 <- NULL

for (i in 1: n_samples) {
  Data_sample <- NULL
  Data_sample <- as.vector(genetic_info[-(1:i),])
  genetic_info_2 <- rbind(genetic_info_2, Data_sample)
}

genetic_info_2 <- as.matrix(genetic_info_2)
dim(genetic_info_2)

#checking dimension of the inputs 
dim(m_input_1) 
dim(m_input_2)
dim(genetic_info_1)
dim(genetic_info_2)

#SEQUENCE BETA VIVERSITY

############################################## MODEL ##############################################


# Parameters --------------------------------------------------------------

input_size <- n_MOTU 
latent_dim <- 2L     
epochs <- 10L            
genetic_dim <- 100L


# Model definition --------------------------------------------------------

# Encoder 1 -----------------------------------------------------------------

layer_genetic_info_1 = layer_input(shape = genetic_dim)

enc_input_1 = layer_input(shape = input_size)
enc_output_1 = enc_input_1 %>% 
  layer_dense(units=150, activation = "relu") %>% 
  layer_dropout(0.2)

h_1 = layer_concatenate(list(enc_output_1, layer_genetic_info_1)) 

output_1_1 = layer_dense(h_1, units=100, activation = "relu") %>% 
  layer_dropout(0.2) %>% 
  layer_dense(units=latent_dim)

output_1_2 = layer_dense(output_1_1, units=50, activation = "relu") %>% 
  layer_dropout(0.2) %>% 
  layer_dense(units=latent_dim)

output_1 = layer_dense(output_1_2, units=20, activation = "relu") %>% 
  layer_dropout(0.2) %>% 
  layer_dense(units=latent_dim)

encoder_1 = keras_model(inputs = c(layer_genetic_info_1, enc_input_1), output_1)
summary(encoder_1)

# Encoder 2 -----------------------------------------------------------------

layer_genetic_info_2 = layer_input(shape = genetic_dim)

enc_input_2 = layer_input(shape = input_size)
enc_output_2 = enc_input_2 %>% 
  layer_dense(units=150, activation = "relu") %>% 
  layer_dropout(0.2) 

h_2 = layer_concatenate(list(enc_output_2 , layer_genetic_info_2)) 

output_2_1 = layer_dense(h_2, units=100, activation = "relu") %>%
  layer_dropout(0.2) %>%
  layer_dense(units=latent_dim)

output_2_2 = layer_dense(output_2_1, units=50, activation = "relu") %>%
  layer_dropout(0.2) %>%
  layer_dense(units=latent_dim)

output_2 = layer_dense(output_2_2, units=20, activation = "relu") %>%
  layer_dropout(0.2) %>%
  layer_dense(units=latent_dim)

encoder_2 = keras_model(inputs = c(layer_genetic_info_2, enc_input_2), output_2)
summary(encoder_2)

# Distances between representation 

euclidian_dist <- function(arg){  
  out_1 <- arg[, 1:(latent_dim)]
  out_2 <- arg[, (latent_dim + 1):(2 * latent_dim)]
  
  
  sqrt(k_sum(k_square(out_1 - out_2), axis = 1))
  
}

eucl_dist = c(enc_output_1, enc_output_2) %>% layer_concatenate() %>% 
  layer_lambda(euclidian_dist)   


Beta_dist_NN <- keras_model(inputs = c(layer_genetic_info_1, enc_input_1, layer_genetic_info_2, enc_input_2), eucl_dist)
summary(Beta_dist_NN)

# Compile and fit ------------------------------------------------------------

#loss function MSE

loss_function<-function(y_pred, y_true){
  k_mean(k_square(y_true - y_pred))
}

Beta_dist_NN %>% compile(optimizer='adam', loss = loss_function)

#Dividing into train and test set

ind <- sample(2:3, nrow(input_1),
              replace = TRUE,
              prob = c(0.8, 0.2))

m_input_1_train <- m_input_1[ind==2,]
m_input_2_train  <- m_input_2[ind==2,]
genetic_info_1_train <- genetic_info_1[ind==2,]
genetic_info_2_train <- genetic_info_2[ind==2,]
beta_values_train <- beta_values[ind==2]

m_input_1_test  <- m_input_1[ind==3,]
m_input_2_test <- m_input_2[ind==3,]
genetic_info_1_test <- genetic_info_1[ind==3,]
genetic_info_2_test <- genetic_info_2[ind==3,]
beta_values_test <- beta_values[ind==3]

# Model training ----------------------------------------------------------

beta_values_train <- as.array(beta_values_train)
dim(beta_values_train)
beta_values_test <- as.array(beta_values_test)
#setting batch size
batch_size = dim(m_input_1_train)[1]

Beta_dist_NN %>% fit(
  x = list(genetic_info_1_train, m_input_1_train, genetic_info_2_train, m_input_2_train), beta_values_train,
  shuffle = FALSE, 
  epochs = epochs, 
  batch_size = dim(m_input_1_train)[1],
  validation_data = list(list(genetic_info_1_test, m_input_1_test, genetic_info_2_test, m_input_2_test), beta_values_test)
)

#predicted values
x_test_encoded_1 <- predict(encoder_1, list( as.matrix(genetic_info), as.matrix(Data[, 1:n_MOTU])), batch_size = 394)
x_test_encoded_2 <- predict(encoder_2, list( as.matrix(genetic_info), as.matrix(Data[, 1:n_MOTU])), batch_size = 394)


# Visualizations -------2D----------------------------------------------

library(ggplot2)
library(dplyr)

motu_richness <- rowSums(Data_presence)


x_test_encoded_1 %>% 
  as_data_frame() %>% 
  ggplot(aes(x = V1, y = V2, col = motu_richness )) + geom_point(aes()) 

x_test_encoded_2 %>% 
  as_data_frame() %>% 
  ggplot(aes(x = V1, y = V2, col = motu_richness )) + geom_point(aes()) 

# Testing -----------------------------------------------------

x_test_encoded <- x_test_encoded_1 
coords_ts <- as.matrix(cbind(x_test_encoded[,1],  x_test_encoded[,2]))
Dist_ls <- (dist(coords_ts, method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
Dist_ls <- Dist_ls / max(Dist_ls)

library(betapart)

jacc_dist <- beta.pair(Data_presence, index.family="jaccard")
matrix_jacc_dist <- as.numeric(jacc_dist$beta.jac)

library("ecodist")

MRM(formula = Dist_ls ~ matrix_jacc_dist, nperm = 1000,
    method = "linear", mrank = FALSE)

#running sequence beta diversity script
MRM(formula = Dist_ls ~ sequence_beta_div , nperm = 1000,
    method = "linear", mrank = FALSE)

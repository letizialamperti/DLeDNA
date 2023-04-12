################################################################################
####                                                                        ####
####          New deep learning-based methods                               ####
####        for visualizing ecosystem properties                            ####
####       using environmental DNA metabarcoding data                       ####
####                                                                        ####
####                                                                        ####
####                                VAESeq                                  ####
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

############################################## MODEL ##############################################


# Parameters --------------------------------------------------------------

original_dim <- n_MOTU 
genetic_dim <- 100L
latent_dim <- 2L   
epochs <- 10   
epsilon_std <- 1.0     

# Model definition --------------------------------------------------------

x_1 <- layer_input(shape = c(original_dim))
x_2 <- layer_input(shape = c(genetic_dim))

h_0 <- layer_dense(object = x_1, units = 150 , activation = "relu") %>%
  layer_dropout(0.2) 
h_1 <- layer_concatenate(list(h_0, x_2))
h <- layer_dense(object = h_1, units = 50, activation = "relu") %>%
  layer_dropout(0.2) 

z_mean <- layer_dense(h, latent_dim)                                   
z_log_var <- layer_dense(h, latent_dim)                                     

sampling <- function(arg){                                                 
  z_mean <- arg[, 1:(latent_dim)]
  z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]
  
  epsilon <- k_random_normal(
    shape = c(k_shape(z_mean)[[1]]), 
    mean=0.,
    stddev=epsilon_std
  )
  
  return(z_mean + k_exp(z_log_var/2)*epsilon)                                   
}

z <- layer_concatenate(c(z_mean, z_log_var)) %>% 
     layer_lambda(sampling)                            

encoder <- keras_model(inputs= c(x_1, x_2), z_mean)
summary(encoder)

#DECODER 
decoder_h <- layer_dense(z, units = 50, activation = "relu")  %>%
  layer_dropout(0.2) 
decoder_h_2 <- layer_dense(decoder_h, units = 150, activation = "relu") %>%
  layer_dropout(0.2) 
decoder_gen_output <- layer_dense(decoder_h_2, units = genetic_dim, activation = "sigmoid", name = 'gen_output') 
decoder_output <- layer_dense(decoder_h_2, units = original_dim, activation = "sigmoid", name = 'main_output')

# end-to-end autoencoder
vae <- keras_model(inputs= c(x_1,x_2), outputs=c(decoder_output,decoder_gen_output))    #modello creato: 1 dense layer 1 che collega il dense al latente, e cosi uguale l'estrazione 

# encoder, from inputs to latent space
encoder <- keras_model(inputs= c(x_1, x_2), z_mean)
summary(encoder)


vae_loss <- function(x, x_decoded_mean){
  xent_loss <- (original_dim/1.0)*loss_binary_crossentropy(x, x_decoded_mean)
  kl_loss <- -0.5*k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L)
  xent_loss + kl_loss
}

vae  %>% compile(
  optimizer = 'adam',
  loss = list(main_output = vae_loss, gen_output = 'binary_crossentropy'),
  loss_weights = list(main_output = 1.0, gen_output = 0.2)
)

#Dividing into train and test set ----------------------------------------------------------

ind <- sample(2:3, nrow(Data),
              replace = TRUE,
              prob = c(0.8, 0.2))

x_train  <- as.matrix(Data[ind==2,])
x_test   <- as.matrix(Data[ind==3,])

genetic_info_train <- as.matrix(genetic_info[ind==2,])
genetic_info_test <- as.matrix(genetic_info[ind==3,])

# Model training ----------------------------------------------------------

vae%>% fit(
  x = list(x_train, genetic_info_train),
  y = list(x_train, genetic_info_train),
  epochs = epochs,
  batch_size = dim(x_train)[1],
  validation_data = list(list(x_test, genetic_info_test), list(x_test, genetic_info_test))
)

#predicted values
x_test_encoded <- predict(encoder, list(as.matrix(Data), as.matrix(genetic_info)), batch_size = batch_size)

# Visualizations -------2D----------------------------------------------

library(ggplot2)
library(dplyr)

Data_presence <- as.matrix((Data>0)+0)
motu_richness <- rowSums(Data_presence)


x_test_encoded %>% 
  as_data_frame() %>% 
  ggplot(aes(x = V1, y = V2, col = motu_richness )) + geom_point(aes())  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "None") 

# Testing -----------------------------------------------------

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

################################################################################
####                                                                        ####
####          New deep learning-based methods                               ####
####        for visualizing ecosystem properties                            ####
####       using environmental DNA metabarcoding data                       ####
####                                                                        ####
####                                                                        ####
####                      Autoencoder genetic info                          ####
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

#import data from http://taraoceans.sb-roscoff.fr/EukDiv/ 
#Database W5. Total V9 rDNA information organized at the OTU level

Database_W5_OTU_occurences <- read.delim("Data/Database_W5_OTU_occurences.tsv")
#selecting taxogroup
Database_W5_OTU_occurences_1 <- Database_W5_OTU_occurences[Database_W5_OTU_occurences$taxogroup == "Dictyochophyceae" , ]

Data <- Database_W5_OTU_occurences_1[, 4:337]

Data <- as.data.frame(t(Data))
#eliminating empty-samples 
Data <- subset(Data, rowSums(Data)>0) 

n_MOTU <- dim(Data)[2]
n_samples <- dim(Data)[1]

Data_presence <- as.matrix((Data>0)+0)

sequences <-  Database_W5_OTU_occurences_1$sequence

### Functions loading
source("Functions_Autoencoder.R")

sequences <-NULL
sequences <- sapply(sequences, uniformize_dimension_seq)
sequences <- as.vector(sequences)

sequences_tran <- NULL
sequences_tran <- lapply(sequences, translate_sequence_fun)
names(sequences_tran) <- sequences

#Input Autoencoder: Tensor of sequences translated and presence/absence 
input_ae <- apply(Data, 1, construct_input_ae)
dim(input_ae)
input_ae <- t(input_ae)

############################################## MODEL ##############################################

# Parameters --------------------------------------------------------------

batch_size <- 100L 
input_size <-  4* max_length_seq * n_MOTU
latent_size <- 100L     
epochs <- 10L            

# Model definition --------------------------------------------------------

# Encoder -----------------------------------------------------------------

enc_input = layer_input(shape = input_size)
enc_output = enc_input %>% 
  layer_dense(units=input_size/32, activation = "relu") %>% 
  layer_dense(units=input_size/64, activation = "relu") %>% 
  layer_dense(units=input_size/128, activation = "relu") %>% 
  layer_dense(units=input_size/256, activation = "relu") %>% 
  layer_dense(units=input_size/512, activation = "relu") %>% 
  layer_dense(units=input_size/1024, activation = "relu") %>% 
  layer_activation_leaky_relu() %>% 
  layer_dense(units=latent_size) %>% 
  layer_activation_leaky_relu()

encoder = keras_model(enc_input, enc_output)
summary(encoder)

# Decoder -----------------------------------------------------------------

dec_input = layer_input(shape = latent_size)
dec_output = dec_input %>% 
  layer_dense(units=input_size/1024, activation = "relu") %>% 
  layer_dense(units=input_size/512, activation = "relu") %>%
  layer_dense(units=input_size/256, activation = "relu") %>% 
  layer_dense(units=input_size/128, activation = "relu") %>%
  layer_dense(units=input_size/64, activation = "relu") %>%
  layer_dense(units=input_size/32, activation = "relu") %>%
  layer_activation_leaky_relu() %>% 
  layer_dense(units = input_size, activation = "sigmoid") %>% 
  layer_activation_leaky_relu()

decoder = keras_model(dec_input, dec_output)

summary(decoder)

# Autoecoder -----------------------------------------------------------------

aen_input = layer_input(shape = input_size)
aen_output = aen_input %>% 
  encoder() %>% 
  decoder()

aen = keras_model(aen_input, aen_output)
summary(aen)

# Compile and fit ------------------------------------------------------------

aen %>% compile(optimizer="adam", loss="binary_crossentropy")

#Dividing into train and test set

#set.seed(111)
ind <- sample(2, nrow(input_ae),
              replace = TRUE,
              prob = c(0.8, 0.2))

x_train  <- input_ae[ind==1,]
x_test   <- input_ae[ind==2,]

# Model training ----------------------------------------------------------

aen %>% fit(
  x_train, x_train, 
  shuffle = TRUE, 
  epochs = epochs, 
  batch_size = batch_size,
  validation_data = list(x_test, x_test)
)

#predicted values
genetic_info <- predict(encoder, input_ae, batch_size = batch_size)
dim(genetic_info)

# Visualizations -------2D----------------------------------------------

library(ggplot2)
library(dplyr)

#PCA on the latent space --> AEgen + PCA 

res.pca <- prcomp((genetic_info), scale = TRUE)

library(betapart)
Data_presence <- as.matrix((Data>0)+0)
motu_richness <- rowSums(Data_presence)

scores <- data.frame(res.pca$x[,1:3])

my_plot <-  scores  %>%
  as_data_frame() %>%
  ggplot(aes(x = PC1 , y = PC2, colour = motu_richness ))  + geom_point() + theme_bw()

my_plot #consider a good plot according to the motu richness gradient

# DLeDNA
New variational autoencoders and deep metric learning based methods improve the dimensionality reduction of environmental DNA metabarcoding data

Here we present two new methods based on deep learning, in which we combine different types of neural networks that we have optimised for eDNA data: the first is based on variational autoencoders (*VAESeq*) and the second on deep metric learning (*ENNBetaDist*). These new methods allow the combination of multiple inputs, in our case the number of sequences found per each molecular operational taxonomic unit (MOTU), with the genetic sequence information of each detected MOTU. 

The genetic autoencoder file contains the part of the model that deals with the latent embedding of sequences. 
The input is the sequences of each MOTU found and the presence/absence of each MOTU in the sample, and the output is the genetic embedding. 
Note that a PCA on the embedding can be performed internally to get an initial visualisation of the data. 

In VAESeq we instead find the part of the model that is based on VAEs. Here we take as input the genetic embedding and the number of sequences until the two-dimensional latent space is generated.

In ENNBetaDist we find the part of the model based on Deep Metric Learning.  We take as input the genetic embedding and the number of sequences until we generate the two-dimensional latent space based on the coupled diversity beta distance.

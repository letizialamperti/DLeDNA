# DLeDNA
ref. doi: https://doi.org/10.1111/1755-0998.13861

New variational autoencoders and deep metric learning based methods improve the dimensionality reduction of environmental DNA metabarcoding data

Here we present two new methods based on deep learning, in which we combine different types of neural networks that we have optimised for eDNA data: the first is based on variational autoencoders (*VAESeq*) and the second on deep metric learning (*ENNBetaDist*). These new methods allow the combination of multiple inputs, in our case the number of sequences found per each molecular operational taxonomic unit (MOTU), with the genetic sequence information of each detected MOTU. 

 ## genetic_autoencoder
 
The genetic autoencoder file contains the part of the model that deals with the latent embedding of sequences. 
The input is the sequences of each MOTU found and the presence/absence of each MOTU in the sample, and the output is the genetic embedding. 
Note that a PCA on the embedding can be performed internally to get an initial visualisation of the data. 

 ## VAESeq 
 
In VAESeq we instead find the part of the model that is based on VAEs. Here we take as input the genetic embedding and the number of sequences until the two-dimensional latent space is generated.

 ## ENNBetaDist 
 
In ENNBetaDist we find the part of the model based on Deep Metric Learning.  We take as input the genetic embedding and the number of sequences until we generate the two-dimensional latent space based on the coupled diversity beta distance.

 ## Data
 
We test our methods on two eukaryotic plankton eDNA datasets from the Tara Ocean expedition (de Vargas et al. 2015). We selected the Dictyochophyceae and Telonemia subsets by taxonomic identification. To help you test the models, we have included the genetic embedding and number of sequences of each subset for easy reference. 

For general use, here is the reference link:
[Database W5. Total V9 rDNA information organized at the OTU level] (http://taraoceans.sb-roscoff.fr/EukDiv/)


Enjoy !

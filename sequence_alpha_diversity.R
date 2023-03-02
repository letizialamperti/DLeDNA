#SEQUENCE ALPHA VIVERSITY

Database_W5_OTU_occurences <- read.delim("Database_W5_OTU_occurences.tsv")
Database_W5_OTU_occurences_1 <- Database_W5_OTU_occurences[Database_W5_OTU_occurences$taxogroup == "Dictyochophyceae" , ]

sequences <- Database_W5_OTU_occurences_1$sequence

max_length_seq <- max(nchar(sequences))

uniformize_dimension_seq <- function(x) (if (nchar(x) < max_length_seq){
  length_temp <- max_length_seq - nchar(x)
  for (i in 1:length_temp) {
    x <- paste(x, 'X', sep="")
  }
  return(x) } else {return(x)})

seq <-NULL
seq <- sapply(sequences, uniformize_dimension_seq)

sequences <- NULL
library(stringr)
for (i in 1:length(seq)) {
  motu_seq <- unlist(str_extract_all(seq[i], boundary("character")))
  sequences <- rbind(sequences, motu_seq)
}

rownames(sequences) <- Database_W5_OTU_occurences_1$sequence

library("ape")

dist_gene <- dist.gene(sequences, method = "percentage", pairwise.deletion = FALSE,
                       variance = FALSE)


library("mFD")
colnames(Data) <- Database_W5_OTU_occurences_1$sequence
dist_gene <- as.matrix(dist_gene)


dim(Data)
dim(dist_gene)

Data_presence <- as.matrix((Data>0)+0)
dim(Data_presence)

seq_alpha_diversity <- alpha.fd.hill( Data_presence,
                                      dist_gene,
                                      q = 0,
                                      tau = "mean",
                                      check_input = TRUE,
                                      details_returned = TRUE)



sequence_alpha_div <- unlist(seq_alpha_diversity[[1]])


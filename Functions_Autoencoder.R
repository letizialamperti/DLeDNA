#############################################################
### Function to genetic Autoencoder



uniformize_dimension_seq <- function(x) (if (nchar(x) < max_length_seq){
  length_temp <- max_length_seq - nchar(x)
  for (i in 1:length_temp) {
    x <- paste(x, 'x', sep="")
  }
  return(x) } else {return(x)})

translate_sequence_fun <- function(x) {
  seq_translated <- NULL 
  for (i in 1:nchar(x)) {
    if(substring(x, i, i) == "a"){
      seq_translated <- rbind(seq_translated, c(1,0,0,0))
    } else if (substring(x, i, i) == "t"){
      seq_translated <- rbind(seq_translated, c(0,1,0,0))
    } else if (substring(x, i, i) == "g"){
      seq_translated <- rbind(seq_translated, c(0,0,1,0))
    } else if (substring(x, i, i) == "c"){
      seq_translated <- rbind(seq_translated, c(0,0,0,1))
    } else if (substring(x, i, i) == "w"){
      seq_translated <- rbind(seq_translated, c(0.5,0,0,0.5))
    } else if (substring(x, i, i) == "s"){
      seq_translated <- rbind(seq_translated, c(0,0.5,0.5,0))
    } else  {
      seq_translated <- rbind(seq_translated, c(0.25,0.25,0.25,0.25))
    }} 
  return(seq_translated)}

construct_input_ae <- function(x) {
  genetic_tensor <- NULL
  for (i in 1:length(x)) {
    if( x[i] != 0){
      current_motu_sequence <- sequences[i]
      matrix_current_motu <-NULL
      matrix_current_motu <- translate_sequence_fun(current_motu_sequence)
      genetic_tensor <- cbind(genetic_tensor, matrix_current_motu)
    }
    else{
      genetic_tensor <- cbind(genetic_tensor, matrix(0, max_length_seq, 4))
    }
  }
  genetic_tensor <- array(genetic_tensor, dim = c(n_MOTU, max_length_seq, 4))
  return(genetic_tensor)
}





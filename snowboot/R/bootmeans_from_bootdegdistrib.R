bootmeans_from_bootdegdistrib <- function(bempd){
  # function takes a dataframe of the bootstrap distribution and returns a
  # vector of sample means for each bootstrap replication. dataframe must have
  # degree values as the column names.
  product_of_k_and_p_k <- matrix(0, dim(bempd)[1],dim(bempd)[2])
  for(i in 1:dim(bempd)[1]){
    product_of_k_and_p_k[i,] <- as.numeric(colnames(bempd))*bempd[i,]
  }
  res <- rowSums(product_of_k_and_p_k)
  res
}

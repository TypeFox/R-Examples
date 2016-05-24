designMat2Vec <-
function(desmat){
  #########
  # Convert from a design matrix (X_t,j=1 is t->j)
  # into a design vector (a map t->j)
  #
  # If a_1...a_m is a vector of effects,
  # then Y_t = desmat_t*a + eps_t and Y_t = desvec(t) + eps_t 
  #########
  desvec = numeric(nrow(desmat))
  for (j in 1:ncol(desmat)){
    desvec[which(desmat[,j]>0)] = j
  }
  return(desvec)
}

# This file contains functions to deal with GAPs.

# The sites of Mu will be prune as GAPs if within cluster data are all GAPs.
prune.Mu <- function(X, X.class, Mu, code.type = .code.type[1]){
  K <- length(unique(X.class)) 
  L <- ncol(X)
  N <- nrow(X)

  if(nrow(Mu) != K){
    stop("X.class and Mu are not consistent.")
  }
  if(ncol(Mu) != L){
    stop("Mu and X are not consistent.")
  }

  if(code.type == "NUCLEOTIDE"){
    GAP <- .nucleotide$nid[.nucleotide$code == "-"]
  } else if(code.type == "SNP"){
    GAP <- .snp$nid[.snp$code == "-"]
  } else{
    stop("The code.type is not implemented.")
  }

  for(k in 1:K){
    sites <- colSums(matrix(X[X.class == k,], ncol = L) == GAP) ==
             sum(X.class == k)
    Mu[k, sites] <- GAP
  }

  Mu
} # End of prune.Mu().


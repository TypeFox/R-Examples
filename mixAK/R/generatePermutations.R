##
##  PURPOSE:   Generate all possible permutations of (1, ..., K)
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   11/02/2010
##
##  FUNCTIONS:  generatePermutations
##
## ======================================================================

## *************************************************************
## generatePermutations
## *************************************************************
generatePermutations <- function(K)
{
  thispackage <- "mixAK"

  if (K <= 0) stop("K must be positive")
  K <- as.integer(K)
  Kfact <- gamma(K + 1)

  RES <- .C("generatePermutations", n_perm    = integer(1),
                                    order     = integer(Kfact * K),
                                    tmp_order = integer(Kfact * K),
                                    rank      = integer(Kfact * K),
                                    K         = as.integer(K),
            PACKAGE = thispackage)

  ##rank <- matrix(RES$rank + 1, ncol=K, byrow=TRUE)
  ##print(rank)
  
  RES <- matrix(RES$order + 1, ncol=K, byrow=TRUE)     ### +1 is here to change C++ indeces starting from 0 to R indeces starting from 1

  ##for (i in 1:Kfact){
  ##  cat("Row ", i, "\n", sep="")
  ##  print(RES[i, rank[i, ]])
  ##  print(rank[i, RES[i, ]])    
  ##}  
  
  return(RES)    
}  


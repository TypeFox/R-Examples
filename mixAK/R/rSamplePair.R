##
##  PURPOSE:   Sample a random pair
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   18/01/2008
##
##  FUNCTIONS:  rSamplePair
##
## ======================================================================

## *************************************************************
## rSamplePair
## *************************************************************
rSamplePair <- function(n, K)
{
  thispackage <- "mixAK"

  if (n <= 0) stop("n must be positive")
  if (K <= 1) stop("K must be higher than 1")

  pairC <- .C("SamplePair_R", j1=integer(n), j2=integer(n), K=as.integer(K), n=as.integer(n), PACKAGE=thispackage)

  pairC$j1 <- pairC$j1 + 1   ### C++ indeces 0,... -> R indeces 1, ...
  pairC$j2 <- pairC$j2 + 1   ### C++ indeces 0,... -> R indeces 1, ...  

  if (n == 1) return(c(pairC$j1, pairC$j2))
  else        return(cbind(pairC$j1, pairC$j2))  
}

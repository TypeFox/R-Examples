##
##  PURPOSE:   Truncated normal distribution
##             * random numbers generation
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   13/11/2007
##
##  FUNCTION:  rTNorm
##
## ======================================================================


## *************************************************************
## rTNorm
## *************************************************************
rTNorm <- function(n, mean=0, sd=1, a, b, trunc)
{
  thispackage <- "mixAK"

  if (n <= 0) stop("n must be positive")
  
  if (length(mean) != length(sd)) stop("mean and sd must have the same length")
  mu.sigma.common <- ifelse(length(mean) == 1, 1, 0)
  if (!mu.sigma.common & length(mean) != n) stop("mean has incorrect length")
  
  if (missing(trunc)){
    trunc <- 4
    a <- b <- 0
    a.b.trunc.common <- 1
  }else{
    a.b.trunc.common <- ifelse(length(trunc) == 1, 1, 0)
    if (a.b.trunc.common){
      if (trunc == 0) b <- 0
      else if (trunc == 1) b <- 0
           else if (trunc == 2) b <- 0
                else if (trunc == 4) a <- b <- 0
    }
    if (any(!(trunc %in% c(0, 1, 2, 3, 4)))) stop("All trunc values must be from {0, 1, 2, 3, 4}.")    
    if (!any(trunc==3)) b <- rep(0, length(trunc))
    if (sum(trunc==4) == length(trunc)) a <- b <- rep(0, length(trunc))    
    if (missing(a) | missing(b)) stop("a and b must be given")
    if (length(trunc) != length(a) | length(trunc) != length(b)) stop("a, b and trunc must have the same length")
    a[trunc==4] <- 0
    b[trunc==0] <- 0
    b[trunc==1] <- 0    
    b[trunc==2] <- 0
    b[trunc==4] <- 0    
    if (any(a[trunc==3] >= b[trunc==3])) stop("a must be lower than b when trunc = 3") 
  }  

  SAMPLE <- .C("rTNorm1_R", x               =double(n),
                            mu              =as.double(mean),
                            sigma           =as.double(sd),
                            a               =as.double(a),
                            b               =as.double(b),
                            trunc           =as.integer(trunc),
                            nx              =as.integer(n),
                            mu.sigma.common =as.integer(mu.sigma.common),
                            a.b.trunc.common=as.integer(a.b.trunc.common),
               PACKAGE=thispackage)

  return(SAMPLE$x)
}  



##
##  PURPOSE:   Calculate initial values of (censored) observations
##             or check supplied inits for consistency.
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    14/02/2010 (by taking sub-code originally included in NMixMCMC function)
##
##  FUNCTIONS:  NMixMCMCinity
##
## ================================================================================================

## *************************************************************
## NMixMCMCinity
## *************************************************************
##
NMixMCMCinity <- function(y0, y1, censor, sd.init, are.Censored, are.Right, are.Exact, are.Left, are.Interval, p, n, inity, random=FALSE)
{
  ##### Calculate initial y if not given  
  if (missing(inity)){
    inity <- y0
    if (are.Censored){
      tmpsd <- matrix(rep(sd.init, n), nrow=n, ncol=p, byrow=TRUE)
      if (random){
        if (are.Right)    inity[censor == 0] <- y0[censor == 0] + abs(rnorm(sum(censor==0), mean=0, sd=tmpsd[censor==0]))
        if (are.Left)     inity[censor == 2] <- y0[censor == 2] - abs(rnorm(sum(censor==2), mean=0, sd=tmpsd[censor==2]))
        if (are.Interval) inity[censor == 3] <- runif(sum(censor==3), min=y0[censor == 3], max=y1[censor == 3])        
      }else{
        if (are.Right)    inity[censor == 0] <- y0[censor == 0] + tmpsd[censor==0]
        if (are.Left)     inity[censor == 2] <- y0[censor == 2] - tmpsd[censor==2]
        if (are.Interval) inity[censor == 3] <- (y0[censor == 3] + y1[censor == 3])/2
      }  
    }  
  }  

  ##### Check and format inity
  if (any(is.na(inity))) stop("NA in init$y")
  if (p == 1) inity <- matrix(inity, ncol=1)
  if (!is.matrix(inity)) stop("init$y must be a matrix")
  if (nrow(inity) != n | ncol(inity) != p) stop("data and init$y mismatch (dimension)")
  if (are.Right) if(any(inity[censor == 0] <= y0[censor == 0])) stop("init$y and y0 mismatch (initial value lower than right-censored observation)")
  if (are.Exact) inity[censor == 1] != y0[censor == 1]
  if (are.Left)  if(any(inity[censor == 2] >= y0[censor == 2])) stop("init$y and y0 mismatch (initial value higher than left-censored observation)")
  if (are.Interval){
    if(any(inity[censor == 3] <= y0[censor == 3])) stop("init$y and y0 mismatch (initial value lower than the left limit of the interval-censored observation)")
    if(any(inity[censor == 3] >= y1[censor == 3])) stop("init$y and y1 mismatch (initial value higher than the right limit of the interval-censored observation)")
  }    

  return(inity)  
}

## code_generic.R --- 
##
## Author: Andreas Kiermeier
##
## Created: 21 Aug 2007 
##
## Purpose: Generic code (mostly) which applies to all types of sampling plans
##          
## Changes:
## 21Aug07: * Created
## ----------------------------------------------------------------------

setGeneric("assess", function(object, PRP, CRP, print=TRUE)
           standardGeneric("assess"))

check.paccept <-
  function(pa){
    ## Purpose: Utility function to check that supplied P(accept) values
    ##          fall within [0,1]
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## pa: a vector of P(accept) values
    ## ----------------------------------------------------------------------
    ## Author: Andreas Kiermeier, Date: 16 May 2007, 10:19
    if (any(pa < 0) | any(pa > 1))
      return(FALSE)
    return(TRUE)
  }

check.quality <-
  function(pd, type){
    ## Purpose: Utility function to check that supplied Proportion defective
    ##          values fall within
    ##          [0,1] for the binomial or hypergeometric
    ##          [0,inf] for the poisson
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## pd: a vector of proportion defective values
    ## ----------------------------------------------------------------------
    ## Author: Andreas Kiermeier, Date: 16 May 2007, 10:19
    if (any(pd < 0))
      return(FALSE)
    if (type %in% c("binomial", "hypergeom", "normal") & any(pd > 1))
      return(FALSE)
    return(TRUE)
  }

## Utility to find k for a given n
find.k <- function(n, pd, pa, interval=c(0,5)){
  tmp <- uniroot(function(x, n, pd, pa){
    pt(x*sqrt(n), df=n-1, ncp=-qnorm(pd)*sqrt(n)) - pa},
                 interval=interval, n=n, pd=pd, pa=1-pa)
  return(tmp$root)
}

find.plan <- function(PRP, CRP,
                      type=c("binomial","hypergeom","poisson","normal"),
                      N,
                      s.type=c("known", "unknown"))
{
  ## Purpose: Find the sampling plan with the smallest sample size, which
  ##          meets a prespecified Producer and Consumer Risk Points.
  ##          
  ##          The convention used here, as in many books, is to use equality
  ##          for the Producer Risk Point rather than the consumer risk point.
  ##
  ##          No consideration is given to "cost functions".
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## PRP   : Producer risk point in the form c(pdefect, paccept)
  ## CRP   : Consumer risk point in the form c(pdefect, paccept)
  ## N     : Population size - only used for hypergeomtric distribution
  ## type  : The distributional assumption
  ## s.type: Only used for 'normal' distribution - indicates whether the
  ##         standard deviation is known or unknown (use sample s.d.)
  ## ----------------------------------------------------------------------
  ## Author: Andreas Kiermeier, Date: 20 Aug 2007, 12:09

  type <- match.arg(type)
  s.type <- match.arg(s.type)

  ## Needs checking that risk points are "valid" - use existing functions
  if (missing(PRP) | missing(CRP))
    stop("Poducer and Consumer Risk Points must be provided.")
  else if(!check.quality(PRP[1], type=type) |
          !check.paccept(PRP[2]) )
    stop("Producer Risk Point - Quality and/or desired P(accept) out of bounds")
  else if(!check.quality(CRP[1], type=type) |
          !check.paccept(CRP[2]) )
    stop("Consumer Risk Point - Quality and/or desired P(accept) out of bounds")
  else if(CRP[1] <= PRP[1])
    stop("Consumer Risk Point quality must be greater than Producer Risk Point quality")

  ## Attributes Sampling Plan - Binomial distribution
  if (type == "binomial") {
    c <- 0
    n <- c+1
    pa1 <- calc.OCbinomial(n=n,c=c,r=c+1,pd=PRP[1])
    pa2 <- calc.OCbinomial(n=n,c=c,r=c+1,pd=CRP[1])

    while (pa2 > CRP[2]) {
      while (pa1 >= PRP[2]) {
        n <- n+1
        pa1 <- calc.OCbinomial(n=n,c=c,r=c+1,pd=PRP[1])
        pa2 <- calc.OCbinomial(n=n,c=c,r=c+1,pd=CRP[1])
        if(pa2 <= CRP[2])
          break
      }
      if(pa2 <= CRP[2] & pa1 >= PRP[2])
          break
      ## No need to reset n - this can stay where it is to speed up
      ## finding a solution.
      c <- c+1
      pa1 <- calc.OCbinomial(n=n,c=c,r=c+1,pd=PRP[1])
      pa2 <- calc.OCbinomial(n=n,c=c,r=c+1,pd=CRP[1])
    }
    return(list(n=n, c=c, r=c+1))
  }
  ## Attributes Sampling Plan - Binomial distribution
  if (type == "hypergeom") {
    c <- 0
    n <- c+1
    pa1 <- calc.OChypergeom(n=n,c=c,r=c+1,N=N,D=PRP[1]*N)
    pa2 <- calc.OChypergeom(n=n,c=c,r=c+1,N=N,D=CRP[1]*N)

    while (pa2 > CRP[2]) {
      while (pa1 >= PRP[2]) {
        n <- n+1
        pa1 <- calc.OChypergeom(n=n,c=c,r=c+1,N=N,D=PRP[1]*N)
        pa2 <- calc.OChypergeom(n=n,c=c,r=c+1,N=N,D=CRP[1]*N)
        if(pa2 <= CRP[2])
          break
      }
      if(pa2 <= CRP[2] & pa1 >= PRP[2])
          break
      ## No need to reset n - this can stay where it is to speed up
      ## finding a solution.
      c <- c+1
      pa1 <- calc.OChypergeom(n=n,c=c,r=c+1,N=N,D=PRP[1]*N)
      pa2 <- calc.OChypergeom(n=n,c=c,r=c+1,N=N,D=CRP[1]*N)
    }
    return(list(n=n, c=c, r=c+1))
  }
  ## Attributes Sampling Plan - Binomial distribution
  if (type == "poisson") {
    c <- 0
    n <- c+1
    pa1 <- calc.OCpoisson(n=n,c=c,r=c+1,pd=PRP[1])
    pa2 <- calc.OCpoisson(n=n,c=c,r=c+1,pd=CRP[1])

    while (pa2 > CRP[2]) {
      while (pa1 >= PRP[2]) {
        n <- n+1
        pa1 <- calc.OCpoisson(n=n,c=c,r=c+1,pd=PRP[1])
        pa2 <- calc.OCpoisson(n=n,c=c,r=c+1,pd=CRP[1])
        if(pa2 <= CRP[2])
          break
      }
      if(pa2 <= CRP[2] & pa1 >= PRP[2])
          break
      ## No need to reset n - this can stay where it is to speed up
      ## finding a solution.
      c <- c+1
      pa1 <- calc.OCpoisson(n=n,c=c,r=c+1,pd=PRP[1])
      pa2 <- calc.OCpoisson(n=n,c=c,r=c+1,pd=CRP[1])
    }
    return(list(n=n, c=c, r=c+1))
  }
  ## Variables Sampling Plan - Normal distribution
  else if (type=="normal") {
    ## With known standard deviation
    if (s.type=="known") {
      n <- ceiling( ((qnorm(1-PRP[2]) + qnorm(CRP[2]))/
                     (qnorm(CRP[1])-qnorm(PRP[1])) )^2)
      k <- qnorm(1-PRP[2])/sqrt(n) - qnorm(PRP[1])
      return(list(n=n, k=k, s.type=s.type))
    }
    ## With unknown standard deviation
    else if (s.type=="unknown") {
      n <- 2 ## Need a minimum of 1 degree of freedom (=n-1) for the NC t-dist
      k <- find.k(n, PRP[1], PRP[2], interval=c(0,1000))
      pa <- 1- pt(k*sqrt(n), df=n-1, ncp=-qnorm(CRP[1])*sqrt(n))
      while(pa > CRP[2]){
        n <- n+1
        k <- find.k(n, PRP[1], PRP[2])
        pa <- 1-pt(k*sqrt(n), df=n-1, ncp=-qnorm(CRP[1])*sqrt(n))
      }
      return(list(n=n, k=k, s.type=s.type))
    }
  }
}



## x1 <- find.plan(c(0.05, 0.95), c(0.15, 0.075), type="bin")
## x <- OC2c(x1$n, x1$c, x1$r, type="bin")
## assess(x, c(0.05, 0.95), c(0.15, 0.075))

## x1 <- find.plan(c(0.05, 0.95), c(0.15, 0.075), type="hyp", N=100)
## x <- OC2c(x1$n, x1$c, x1$r, type="hyp", N=100)
## assess(x, c(0.05, 0.95), c(0.15, 0.075))

## x1 <- find.plan(c(0.05, 0.95), c(0.15, 0.075), type="pois")
## x <- OC2c(x1$n, x1$c, x1$r, type="pois")
## assess(x, c(0.05, 0.95), c(0.15, 0.075))


## The following examples come from Guenther's book

## PRP <- c(0.01, 0.95)
## CRP <- c(0.10, 0.1)
## x1 <- find.plan(PRP=PRP, CRP=CRP, type="nor", s.type="unknown")
## x <- OCvar(x1$n, x1$k, s.type=x1$s.type, pd=seq(0,0.2, by=0.002))
## plot(x)
## points(PRP[1], PRP[2], col="red", pch=19); points(CRP[1], CRP[2], col="red", pch=19)
## assess(x, PRP=PRP, CRP=CRP)

## PRP <- c(0.05, 0.95)
## CRP <- c(0.20, 0.1)
## x1 <- find.plan(PRP=PRP, CRP=CRP, type="nor", s.type="known")
## x <- OCvar(x1$n, x1$k, s.type=x1$s.type, pd=seq(0,0.2, by=0.002))
## plot(x)
## points(PRP[1], PRP[2], col="red", pch=19); points(CRP[1], CRP[2], col="red", pch=19)
## assess(x, PRP=PRP, CRP=CRP)





### Local Variables:
### comment-start: "## "
### fill-column: 80
### End:


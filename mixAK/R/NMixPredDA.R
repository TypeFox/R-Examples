##
##  PURPOSE:    Discriminant analysis based on mixture model
##              Determine component probabilities for "new" observations
##              based on sampled mixture components
##              (based on NMixMCMC)
##              * "new" observations are allowed to be censored
##              * sample is re-labeled using the 'order' component before computation of allocation probabilities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   13/02/2010 
##
##  FUNCTION:  NMixPredDA
##             
##
## ======================================================================

## *************************************************************
## NMixPredDA
## *************************************************************
NMixPredDA <- function(object, y0, y1, censor, inity, info)
{
  thispackage <- "mixAK"
  
  ########## Check MCMC object, parameters of MCMC
  ########## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (class(object) != "NMixMCMC") stop("object must be of class NMixMCMC")
  if (object$prior$priorK != "fixed") stop("number of mixture components was not fixed")

  if (object$nx_w > 1) stop("This function has not (yet) been implemented if a factor covariate on mixture weights is present.")
  
  if (is.null(object$K) | is.null(object$w) | is.null(object$mu) | is.null(object$Li) | is.null(object$Q) | is.null(object$Sigma)){
    stop("object does not contain sampled values")
  }  

  keepMCMC <- length(object$w) / object$K[1]
  if (missing(info)) info <- keepMCMC
  if (info <= 0 | info > keepMCMC) info <- keepMCMC
  

  ########## Data if not given (taken from object)
  ########## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (missing(y0)){
    if (object$dim == 1){
      y0 <- object$scale$scale * object$Cpar$z0 + object$scale$shift
      y1 <- object$scale$scale * object$Cpar$z1 + object$scale$shift
    }else{
      y0 <- matrix(rep(object$scale$scale, object$Cpar$dimy["n"]), ncol=object$dim, byrow=TRUE) * object$Cpar$z0 + matrix(rep(object$scale$shift, object$Cpar$dimy["n"]), ncol=object$dim, byrow=TRUE)
      y1 <- matrix(rep(object$scale$scale, object$Cpar$dimy["n"]), ncol=object$dim, byrow=TRUE) * object$Cpar$z1 + matrix(rep(object$scale$shift, object$Cpar$dimy["n"]), ncol=object$dim, byrow=TRUE)
    }
    censor <- object$Cpar$censor

    if (missing(inity)){
      inity <- object$state.first$y
    }
  }

  
  ########## Data
  ########## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  dd <- NMixMCMCdata(y0=y0, y1=y1, censor=censor)
  rm(list=c("y0", "y1", "censor"))
  if (dd$p != object$dim) stop("supplied y0 has another dimension than the data used to generate MCMC sample")
  z0    <- (dd$y0 - matrix(rep(object$scale$shift, dd$n), ncol=object$dim, byrow=TRUE)) / matrix(rep(object$scale$scale, dd$n), ncol=object$dim, byrow=TRUE)
  z1    <- (dd$y1 - matrix(rep(object$scale$shift, dd$n), ncol=object$dim, byrow=TRUE)) / matrix(rep(object$scale$scale, dd$n), ncol=object$dim, byrow=TRUE)
  
  
  ########## Temporar initial values (will be used to generate inity if not given)
  ########## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  tmpinity <- dd$y0
  if (dd$are.Interval) tmpinity[dd$censor == 3] <- (dd$y0[dd$censor == 3] + dd$y1[dd$censor == 3])/2      


  ########## Initial y and z
  ########## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (missing(inity)){
    inity <- NMixMCMCinity(y0=dd$y0, y1=dd$y1, censor=dd$censor, sd.init=sd(tmpinity),
                           are.Censored=dd$are.Censored, are.Right=dd$are.Right, are.Exact=dd$are.Exact, are.Left=dd$are.Left, are.Interval=dd$are.Interval,
                           p=dd$p, n=dd$n, random=FALSE)
  }else{
    inity <- NMixMCMCinity(y0=dd$y0, y1=dd$y1, censor=dd$censor, sd.init=sd(tmpinity),
                           are.Censored=dd$are.Censored, are.Right=dd$are.Right, are.Exact=dd$are.Exact, are.Left=dd$are.Left, are.Interval=dd$are.Interval,
                           p=dd$p, n=dd$n, inity=inity)
  }
  initz <- (inity - matrix(rep(object$scale$shift, dd$n), ncol=object$dim, byrow=TRUE)) / matrix(rep(object$scale$scale, dd$n), ncol=object$dim, byrow=TRUE)
  

  ########## Main computation
  ########## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MCMC <- .C("NMix_PredDA",
             y0           = as.double(t(z0)),
             y1           = as.double(t(z1)),
             censor       = as.integer(t(dd$censor)),
             dimy         = as.integer(c(dd$p, dd$n)),
             keepMCMC     = as.integer(keepMCMC),
             info         = as.integer(info),
             K            = as.integer(object$K[1]),
             chw          = as.double(t(object$w)),
             chmu         = as.double(t(object$mu)),
             chSigma      = as.double(t(object$Sigma)),
             chLi         = as.double(t(object$Li)),
             chrank       = as.integer(t(object$rank - 1)),
             y            = as.double(t(initz)),
             r            = integer(dd$n),
             sum_Ir       = integer(dd$n * object$K[1]),
             hatPr_y      = double(dd$n * object$K[1]),
             err          = as.integer(0),
             PACKAGE = thispackage)
  if (MCMC$err) stop("Something went wrong.")
  

  ########## Discriminant analysis based on posterior P(alloc = k | y) or on P(alloc = k | theta, b, y) 
  ########## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (object$K[1] == 1){
    poster.comp.prob_u <- poster.comp.prob_b <- matrix(1, nrow = dd$n, ncol = 1)
  }else{

    ### Using mean(I(r=k))
    MCMC$sum_Ir <- matrix(MCMC$sum_Ir, ncol = object$K[1], nrow = dd$n, byrow = TRUE)
    Denom <- apply(MCMC$sum_Ir, 1, sum)       ### this should be a vector of length n with all elements equal to the number of saved MCMC iterations 
    poster.comp.prob_u <- MCMC$sum_Ir / matrix(rep(Denom, object$K[1]), ncol = object$K[1], nrow = dd$n)

    ### Using mean(P(r=k | theta, b, y))
    poster.comp.prob_b <- matrix(MCMC$hatPr_y, ncol = object$K[1], nrow = dd$n, byrow = TRUE)
  }  


  ########## Use poster.comp.prob_b to discriminate and return
  ########## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  colnames(poster.comp.prob_b) <- paste("prob", 1:object$K[1], sep="")  
  RET <- as.data.frame(poster.comp.prob_b)
  RET[,"component"] <- apply(poster.comp.prob_b, 1, which.max)  
  return(RET)
}
  

##
##  PURPOSE:   Handle prior.b argument of GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  LOG:        20090805
##
##  FUNCTIONS:  GLMM_MCMCprior.b
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCprior.b
## *************************************************************
##
GLMM_MCMCprior.b <- function(prior.b, scale.b, dimb, iEranefVec, iSDranefVec)
{
##### Variables in the resulting object:
#####             CpriorInt_b
#####             CpriorDouble_b
##### -----------------------------------------------------------------------------------------------------------          
  if (dimb){
    LTb <- (dimb * (dimb + 1))/2
    Imat <- diag(dimb)
    rowsI <- row(Imat)[lower.tri(row(Imat), diag=TRUE)]
    colsI <- col(Imat)[lower.tri(col(Imat), diag=TRUE)] 
    naamLTb <- paste(".", rowsI, ".", colsI, sep="")    
    
    if (missing(iEranefVec))  iEranefVec <- rep(0, dimb)
    if (missing(iSDranefVec)) iSDranefVec <- rep(1, dimb)    
    
    iEbscaled <- (iEranefVec - scale.b$shift) / scale.b$scale
    iSDbscaled <- iSDranefVec / scale.b$scale    
    Rbb <- 6*iSDbscaled
    bbmin <- iEbscaled - 3*iSDbscaled
    bbmax <- iEbscaled + 3*iSDbscaled    
    if (dimb == 1) bbVar <- iSDbscaled^2
    else           bbVar <- diag(iSDbscaled^2)
    
    if (missing(prior.b)) prior.b <- list()
    if (!is.list(prior.b)) stop("prior.b must be a list")

    inprior.b <- names(prior.b)
    ib.distribution <- match("distribution", inprior.b, nomatch=NA)    
    ib.priorK       <- match("priorK",       inprior.b, nomatch=NA)
    ib.priormuQ     <- match("priormuQ",     inprior.b, nomatch=NA)
    ib.Kmax         <- match("Kmax",         inprior.b, nomatch=NA)
    ib.lambda       <- match("lambda",       inprior.b, nomatch=NA)
    ib.delta        <- match("delta",        inprior.b, nomatch=NA)
    ib.xi           <- match("xi",           inprior.b, nomatch=NA)
    ib.ce           <- match("ce",           inprior.b, nomatch=NA)
    ib.D            <- match("D",            inprior.b, nomatch=NA)
    ib.zeta         <- match("zeta",         inprior.b, nomatch=NA)
    ib.gD           <- match("gD",           inprior.b, nomatch=NA)
    ib.hD           <- match("hD",           inprior.b, nomatch=NA)
    ib.gdf          <- match("gdf",          inprior.b, nomatch=NA)
    ib.hdf          <- match("hdf",          inprior.b, nomatch=NA)                
    

    ##### prior.b:  distribution
    ##### -----------------------------------------------
    if (is.na(ib.distribution)) prior.b$distribution <- "normal"
    if (length(prior.b$distribution) != 1) stop("prior.b$distribution must be of length 1")
    Cbdistribution <- pmatch(prior.b$distribution, table=c("normal", "MVT"), nomatch=0) - 1
    if (Cbdistribution == -1) stop("prior.b$distribution must be one of normal/MVT")
    
    ##### prior.b:  priorK
    ##### -----------------------------------------------
    if (is.na(ib.priorK)) prior.b$priorK <- "fixed"
    if (length(prior.b$priorK) != 1) stop("prior.b$priorK must be of length 1")
    CbpriorK <- pmatch(prior.b$priorK, table=c("fixed", "uniform", "tpoisson"), nomatch=0) - 1
    if (CbpriorK == -1) stop("prior.b$priorK must be one of fixed/uniform/tpoisson")
    if (CbpriorK > 0) stop("prior.b$priorK other than fixed not (yet) implemented")

    ##### prior.b:  priormuQ
    ##### -----------------------------------------------
    if (is.na(ib.priormuQ)) prior.b$priormuQ <- "independentC"
    if (length(prior.b$priormuQ) != 1) stop("prior.b$priormuQ must be of length 1")  
    CbpriormuQ <- pmatch(prior.b$priormuQ, table=c("naturalC", "independentC"), nomatch=0) - 1
    if (CbpriormuQ == -1) stop("prior.b$priormuQ must be one of naturalC/independentC")
  
    ##### prior.b:  Kmax
    ##### -----------------------------------------------
    if (is.na(ib.Kmax)) prior.b$Kmax <- 1
    if (length(prior.b$Kmax) != 1) stop("prior.b$Kmax must be of length 1")
    if (is.na(prior.b$Kmax)) stop("NA in prior.b$Kmax")
    if (prior.b$Kmax <= 0) stop("prior.b$Kmax must be positive")
    CbKmax <- as.numeric(prior.b$Kmax)
    
    ##### prior.b:  lambda
    ##### -----------------------------------------------
    if (CbpriorK == 2){    ## truncated Poisson prior for K
      if (is.na(ib.lambda)) stop("prior.b$lambda must be given when prior.b$priorK = tpoisson")
      if (length(prior.b$lambda) != 1) stop("prior.b$lambda must be of length 1")
      if (is.na(prior.b$lambda)) stop("NA in prior.b$lambda")
      if (prior.b$lambda <= 0) stop("prior.b$lambda must be positive")    
    }else{
      prior.b$lambda <- 0
    }
    Cblambda <- as.numeric(prior.b$lambda)
    names(Cblambda) <- "lambda"

    ##### prior.b:  delta
    ##### -----------------------------------------------
    if (is.na(ib.delta)) prior.b$delta <- 1
    if (length(prior.b$delta) != 1) stop("prior.b$delta must be of length 1")
    if (is.na(prior.b$delta)) stop("NA in prior.b$delta")  
    if (prior.b$delta <= 0) stop("prior.b$delta must be positive")
    Cbdelta <- as.numeric(prior.b$delta)
    names(Cbdelta) <- "delta"

    ##### prior.b:  xi
    ##### -----------------------------------------------
    if (is.na(ib.xi)) prior.b$xi <- matrix(rep(iEbscaled, CbKmax), nrow=CbKmax, ncol=dimb, byrow=TRUE)
    if (any(is.na(prior.b$xi))) stop("NA in prior.b$xi")  
    if (dimb == 1){
      if (length(prior.b$xi) == 1) prior.b$xi <- rep(prior.b$xi, CbKmax)                                                       ## common prior.b$xi for all mixture components
      if (length(prior.b$xi) != CbKmax) stop(paste("prior.b$xi must be of length ", CbKmax, sep=""))
      prior.b$xi <- as.numeric(prior.b$xi)
      names(prior.b$xi) <- paste("xi", 1:CbKmax, sep="")    
      Cbxi <- prior.b$xi
    }else{
      if (length(prior.b$xi) == dimb) prior.b$xi <- matrix(rep(as.numeric(prior.b$xi), each=CbKmax), nrow=CbKmax, ncol=dimb)   ## common prior.b$xi for all mixture components
      if (CbKmax == 1) prior.b$xi <- matrix(as.numeric(prior.b$xi), nrow=1)    
      if (!is.matrix(prior.b$xi)) stop("prior.b$xi must be a matrix")
      if (ncol(prior.b$xi) != dimb) stop(paste("prior.b$xi must have ", dimb, " columns", sep=""))
      if (nrow(prior.b$xi) != CbKmax) stop(paste("prior.b$xi must have ", CbKmax, " rows", sep=""))
      rownames(prior.b$xi) <- paste("j", 1:CbKmax, sep="")
      colnames(prior.b$xi) <- paste("m", 1:dimb, sep="")
      Cbxi <- as.numeric(t(prior.b$xi))
      names(Cbxi) <- paste("xi", rep(1:CbKmax, each=dimb), ".", rep(1:dimb, CbKmax), sep="")
    }  
    if (any(is.na(Cbxi))) stop("NA in prior.b$xi")    

    ##### prior.b:  ce
    ##### -----------------------------------------------
    if (CbpriormuQ == 0){    ## natural conjugate prior for (mu, Q)
      if (is.na(ib.ce)) prior.b$ce <- rep(1, CbKmax)
      if (length(prior.b$ce) == 1) prior.b$ce <- rep(prior.b$ce, CbKmax)
      if (length(prior.b$ce) != CbKmax) stop(paste("prior.b$ce must be of length ", CbKmax, sep=""))
      if (any(is.na(prior.b$ce))) stop("NA in prior.b$ce")
      if (any(prior.b$ce <= 0)) stop("prior.b$ce must be positive")
      prior.b$ce <- as.numeric(prior.b$ce)
    }else{
      prior.b$ce <- rep(0, CbKmax)
    }  
    Cbce <- prior.b$ce
    names(Cbce) <- names(prior.b$ce) <- paste("c", 1:CbKmax, sep="")

    ##### prior.b:  D
    ##### -----------------------------------------------
    if (CbpriormuQ == 1){    ## independent conjugate prior for (mu, Q)
      if (is.na(ib.D)){
        if (dimb == 1) prior.b$D <- rep(Rbb^2, CbKmax)
        else           prior.b$D <- t(matrix(rep(diag(Rbb^2), CbKmax), nrow=dimb, ncol=CbKmax*dimb))
      }
      if (any(is.na(prior.b$D))) stop("NA in prior.b$D")    
      if (dimb == 1){
        if (length(prior.b$D) == 1) prior.b$D <- rep(prior.b$D, CbKmax)
        if (length(prior.b$D) != CbKmax) stop(paste("prior.b$D must be of length ", CbKmax, sep=""))
        prior.b$D <- as.numeric(prior.b$D)
        names(prior.b$D) <- paste("D", 1:CbKmax, sep="")      
        if (any(prior.b$D <= 0)) stop("prior.b$D must be positive")
        CbDinv <- 1/prior.b$D
        names(CbDinv) <- paste("Dinv", 1:CbKmax, sep="")      
      }else{
        if (!is.matrix(prior.b$D)) stop("prior.b$D must be a matrix")
        if (ncol(prior.b$D) != dimb) stop(paste("prior.b$D must have ", dimb, " columns", sep=""))
        if (nrow(prior.b$D) == dimb){
          if (any(prior.b$D[lower.tri(prior.b$D)] != t(prior.b$D)[lower.tri(prior.b$D)])) stop("prior.b$D must be a symmetric matrix")
          err <- try(Dinv <- chol(prior.b$D), silent=TRUE)
          if (class(err) == "try-error") stop("Cholesky decomposition of prior.b$D failed")
          Dinv <- chol2inv(Dinv)
          CbDinv <- rep(Dinv[lower.tri(Dinv, diag=TRUE)], CbKmax)        
          prior.b$D <- matrix(rep(as.numeric(t(prior.b$D)), CbKmax), nrow=dimb*CbKmax, ncol=dimb, byrow=TRUE)
        }else{  
          if (nrow(prior.b$D) != CbKmax*dimb) stop(paste("prior.b$D must have ", CbKmax, " times ", dimb, " rows", sep=""))
          CbDinv <- numeric(0)
          for (j in 1:CbKmax){
            Dinv <- prior.b$D[((j-1)*dimb+1):(j*dimb),]
            if (any(Dinv[lower.tri(Dinv)] != t(Dinv)[lower.tri(Dinv)])) stop(paste(j, "-th block of prior.b$D is not symmetric", sep=""))
            err <- try(Dinv <- chol(Dinv), silent=TRUE)
            if (class(err) == "try-error") stop(paste("Cholesky decomposition of the ", j, "-th block of prior.b$D failed", sep=""))
            Dinv <- chol2inv(Dinv)
            CbDinv <- c(CbDinv, Dinv[lower.tri(Dinv, diag=TRUE)])
          }  
        }
        colnames(prior.b$D) <- paste("m", 1:dimb, sep="")
        rownames(prior.b$D) <- paste("j", rep(1:CbKmax, each=dimb), ".", rep(1:dimb, CbKmax), sep="")
        names(CbDinv) <- paste("Dinv", rep(1:CbKmax, each=LTb), rep(naamLTb, CbKmax), sep="")
      }  
    }else{
      if (dimb == 1){
        prior.b$D <- rep(1, CbKmax)
        names(prior.b$D) <- paste("D", 1:CbKmax, sep="")
        CbDinv <- 1/prior.b$D
        names(CbDinv) <- paste("Dinv", 1:CbKmax, sep="")      
      }else{
        prior.b$D <- matrix(rep(as.numeric(diag(dimb)), CbKmax), nrow=dimb*CbKmax, ncol=dimb, byrow=TRUE)
        colnames(prior.b$D) <- paste("m", 1:dimb, sep="")
        rownames(prior.b$D) <- paste("j", rep(1:CbKmax, each=dimb), ".", rep(1:dimb, CbKmax), sep="")

        Dinv <- diag(dimb)
        CbDinv <- rep(Dinv[lower.tri(Dinv, diag=TRUE)], CbKmax)
        names(CbDinv) <- paste("Dinv", rep(1:CbKmax, each=LTb), rep(naamLTb, CbKmax), sep="")      
      }  
    }  
        
    ##### prior.b:  zeta
    ##### -----------------------------------------------
    if (is.na(ib.zeta)) prior.b$zeta <- dimb + 1
    if (length(prior.b$zeta) != 1) stop("prior.b$zeta must be of length 1")  
    if (is.na(prior.b$zeta)) stop("NA in prior.b$zeta")
    if (prior.b$zeta <= dimb - 1) stop(paste("prior.b$zeta must be higher than ", dimb - 1, sep=""))
    Cbzeta <- as.numeric(prior.b$zeta)
    names(Cbzeta) <- "zeta"

    ##### prior.b:  gD
    ##### -----------------------------------------------
    if (is.na(ib.gD)) prior.b$gD <- rep(0.2, dimb)
    if (length(prior.b$gD) == 1) prior.b$gD <- rep(prior.b$gD, dimb)
    if (length(prior.b$gD) != dimb) stop(paste("prior.b$gD must be of length ", dimb, sep=""))  
    if (any(is.na(prior.b$gD))) stop("NA in prior.b$gD")
    if (any(prior.b$gD <= 0)) stop("prior.b$gD must be positive")
    CbgD <- as.numeric(prior.b$gD)
    names(CbgD) <- paste("g", 1:dimb, sep="")

    ##### prior.b:  hD
    ##### -----------------------------------------------
    if (is.na(ib.hD)) prior.b$hD <- 10/(Rbb^2)
    if (length(prior.b$hD) == 1) prior.b$hD <- rep(prior.b$hD, dimb)
    if (length(prior.b$hD) != dimb) stop(paste("prior.b$hD must be of length ", dimb, sep=""))  
    if (any(is.na(prior.b$hD))) stop("NA in prior.b$hD")
    if (any(prior.b$hD <= 0)) stop("prior.b$hD must be positive")
    CbhD <- as.numeric(prior.b$hD)
    names(CbhD) <- paste("h", 1:dimb, sep="")

    ##### prior.b:  gdf
    ##### -----------------------------------------------
    if (is.na(ib.gdf)) prior.b$gdf <- 1
    if (length(prior.b$gdf) != 1) prior.b$gdf <- prior.b$gdf[1]
    if (any(is.na(prior.b$gdf))) stop("NA in prior.b$gdf")
    if (any(prior.b$gdf <= 0)) stop("prior.b$gdf must be positive")
    Cbgdf <- as.numeric(prior.b$gdf)
    names(Cbgdf) <- "gdf"

    ##### prior.b:  hdf
    ##### -----------------------------------------------
    if (is.na(ib.hdf)) prior.b$hdf <- 0.005
    if (length(prior.b$hdf) != 1) prior.b$hdf <- prior.b$hdf[1]
    if (any(is.na(prior.b$hdf))) stop("NA in prior.b$hdf")
    if (any(prior.b$hdf <= 0)) stop("prior.b$hdf must be positive")
    Cbhdf <- as.numeric(prior.b$hdf)
    names(Cbhdf) <- "hdf"
    
    ##### put all together
    ##### -----------------------------------------------
    CpriorInt_b <- c(Cbdistribution, CbpriorK, CbpriormuQ, CbKmax)
    names(CpriorInt_b) <- c("distribution", "priorK", "priormuQ", "Kmax")  
    CpriorDouble_b<- c(Cblambda, Cbdelta, Cbxi, Cbce, CbDinv, Cbzeta, CbgD, CbhD, Cbgdf, Cbhdf)
  }else{
    prior.b <- list(priorK="fixed", priormuQ="independentC", Kmax=0,
                    lambda=0, delta=0, xi=0, ce=0, D=0, zeta=0, gD=0, hD=0, gdf=0, hdf=0)
    CpriorInt_b <- c(0, 1, 0)
    CpriorDouble_b <- rep(0, 10)
  }  

  RET <- list(prior.b        = prior.b,
              CpriorInt_b    = CpriorInt_b,
              CpriorDouble_b = CpriorDouble_b)
  return(RET)
}


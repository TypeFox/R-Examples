##
##  PURPOSE:   Handle init.b or init2.b argument of GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    02/11/2011    
##
##  FUNCTIONS:  GLMM_MCMCinit.b
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCinit.b
## *************************************************************
##
GLMM_MCMCinit.b <- function(init.b, prior.b, scale.b, id, dimb, LTb, naamLTb, I, ibMat, iEranefVec, iSDranefVec, number="")
{
  if (dimb){
    if (missing(init.b)) init.b <- list()
    if (!is.list(init.b)) stop("init.b must be a list")

    ininit.b    <- names(init.b)
    ib.b        <- match("b",        ininit.b, nomatch=NA)
    ib.K        <- match("K",        ininit.b, nomatch=NA)
    ib.w        <- match("w",        ininit.b, nomatch=NA)
    ib.mu       <- match("mu",       ininit.b, nomatch=NA)
    ib.Sigma    <- match("Sigma",    ininit.b, nomatch=NA)
    ib.Li       <- match("Li",       ininit.b, nomatch=NA)
    ib.gammaInv <- match("gammaInv", ininit.b, nomatch=NA)
    ib.df       <- match("df",       ininit.b, nomatch=NA)    
    ib.r        <- match("r",        ininit.b, nomatch=NA)

    ##### init.b:  b (not scaled and not shifted!!!)
    ##### ----------------------------------------------------  
    if (is.na(ib.b)) init.b$b <- as.matrix(ibMat)
    if (dimb == 1) init.b$b <- matrix(as.numeric(init.b$b), ncol=1)
    if (is.data.frame(init.b$b)) init.b$b <- as.matrix(init.b$b)
    if (!is.matrix(init.b$b)) stop(paste("init", number, ".b$b must be a matrix", sep=""))
    if (ncol(init.b$b) != dimb) stop(paste("init", number, ".b$b must have ", dimb, " columns", sep=""))
    if (nrow(init.b$b) != I) stop(paste("init", number, ".b$b must have ", I, " rows", sep=""))
    if (is.null(rownames(init.b$b))) rownames(init.b$b) <- unique(id)
    if (is.null(colnames(init.b$b))) colnames(init.b$b) <- paste("b", 1:dimb, sep="")
    if (any(is.na(init.b$b))) stop(paste("NA in init", number, ".b$b", sep=""))
    
    ##### init.b:  K
    ##### ----------------------------------------------------  
    if (is.na(ib.K)){
      if (prior.b$priorK == "fixed") init.b$K <- prior.b$Kmax
      else                           init.b$K <- 1
    }
    if (prior.b$priorK == "fixed") init.b$K <- prior.b$Kmax
    if (length(init.b$K) != 1) stop(paste("init", number, ".b$K must be of length 1", sep=""))
    if (is.na(init.b$K)) stop(paste("NA in init", number, ".b$K", sep=""))
    if (init.b$K <= 0 | init.b$K > prior.b$Kmax) stop(paste("init", number, ".b$K out of the range", sep=""))
    
    ##### init.b:  w
    ##### ----------------------------------------------------  
    if (is.na(ib.w)){
      init.b$w <- rep(1, init.b$K)/init.b$K
    }  
    init.b$w <- as.numeric(init.b$w)
    if (length(init.b$w) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$w <- init.b$w[1:init.b$K]  
    names(init.b$w) <- paste("w", 1:init.b$K, sep="")
    if (any(is.na(init.b$w))) stop(paste("NA in init", number, ".b$w", sep=""))  
    if (length(init.b$w) != init.b$K) stop(paste("init", number, ".b$w must be of length ", init.b$K, sep=""))
    if (any(init.b$w < 0)) stop(paste("init", number, ".b$w may not be negative", sep=""))
    init.b$w <- init.b$w / sum(init.b$w)

    ##### init.b:  mu
    ##### ----------------------------------------------------
    Rbb <- 6 * (iSDranefVec / scale.b$scale)
    bbmin <- (iEranefVec - scale.b$shift) / scale.b$scale - 0.5 * Rbb
    bbmax <- (iEranefVec - scale.b$shift) / scale.b$scale + 0.5 * Rbb
    
    if (is.na(ib.mu)){      
      if (dimb == 1){
        afstand <- Rbb/(init.b$K + 1)
        init.b$mu <- seq(bbmin+afstand, bbmax-afstand, length=init.b$K)
      }else{
        afstand <- Rbb/(init.b$K + 1)
        init.b$mu <- matrix(NA, nrow=init.b$K, ncol=dimb)
        for (j in 1:dimb) init.b$mu[,j] <- seq(bbmin[j]+afstand[j], bbmax[j]-afstand[j], length=init.b$K)
      }  
    }
    if (any(is.na(init.b$mu))) stop(paste("NA in init", number, ".b$mu", sep=""))          
    if (dimb == 1){
      init.b$mu <- as.numeric(init.b$mu)
      if (length(init.b$mu) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$mu <- init.b$mu[1:init.b$K]          
      if (length(init.b$mu) != init.b$K) stop(paste("init", number, ".b$mu must be of length ", init.b$K, sep=""))
      names(init.b$mu) <- paste("mu", 1:init.b$K, sep="")    
    }else{
      if (!is.matrix(init.b$mu)) stop(paste("init", number, ".b$mu must be a matrix", sep=""))
      if (ncol(init.b$mu) != dimb) stop(paste("init", number, ".b$mu must have ", dimb, " columns", sep=""))
      if (nrow(init.b$mu) != init.b$K) stop(paste("init", number, ".b$mu must have ", init.b$K, " rows", sep=""))
      rownames(init.b$mu) <- paste("j", 1:init.b$K, sep="")
      colnames(init.b$mu) <- paste("m", 1:dimb, sep="")        
    }

    ##### init.b:  Sigma and Li
    ##### ----------------------------------------------------
    if (dimb == 1) bbVar <- (iSDranefVec / scale.b$scale)^2
    else           bbVar <- diag((iSDranefVec / scale.b$scale)^2)
    
    if (is.na(ib.Sigma)){            
      if (is.na(ib.Li)){       ### Sigma and Li are computed from the data
        if (dimb == 1){
          init.b$Sigma <- rep(bbVar, init.b$K)
          names(init.b$Sigma) <- paste("Sigma", 1:init.b$K, sep="")
          init.b$Li <- sqrt(1 / init.b$Sigma)
          names(init.b$Li) <- paste("Li", 1:init.b$K, sep="")      
        }else{
          init.b$Sigma <- matrix(rep(t(bbVar), init.b$K), ncol=dimb, byrow=TRUE)
          Sigmainv <- chol(bbVar)        
          Sigmainv <- chol2inv(Sigmainv)
          Litmp <- t(chol(Sigmainv))
          init.b$Li <- rep(Litmp[lower.tri(Litmp, diag=TRUE)], init.b$K)
          rownames(init.b$Sigma) <- paste("j", rep(1:init.b$K, each=dimb), ".", rep(1:dimb, init.b$K), sep="")
          colnames(init.b$Sigma) <- paste("m", 1:dimb, sep="")                
          names(init.b$Li) <- paste("Li", rep(1:init.b$K, each=LTb), rep(naamLTb, init.b$K), sep="")
        }        
      }else{                 ### Li is checked and Sigma is computed from Li
        if (any(is.na(init.b$Li))) stop(paste("NA in init", number, ".b$Li", sep="")) 
        if (dimb == 1){
          if (length(init.b$Li) == 1) init.b$Li <- rep(init.b$Li, init.b$K)
          if (length(init.b$Li) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$Li <- init.b$Li[1:init.b$K]
          if (length(init.b$Sigma) != init.b$K) stop(paste("init", number, ".b$Sigma must be of length ", init.b$K, sep=""))
          init.b$Li <- as.numeric(init.b$Li)
          names(init.b$Li) <- paste("Li", 1:init.b$K, sep="")      
          if (any(init.b$Li <= 0)) stop(paste("init", number, ".b$Li must be positive", sep=""))
          init.b$Sigma <- (1 / init.b$Li)^2
          names(init.b$Sigma) <- paste("Sigma", 1:init.b$K, sep="")      
        }else{
          if (length(init.b$Li) == LTb){
            tmpSigma <- matrix(0, nrow=dimb, ncol=dimb)
            tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- init.b$Li
            tmpSigma <- tmpSigma %*% t(tmpSigma)
            err <- try(tmpSigma <- chol(tmpSigma), silent=TRUE)
            if (class(err) == "try-error") stop(paste("init", number, ".b$Li does not lead to a positive definite matrix", sep=""))
            tmpSigma <- chol2inv(tmpSigma)
            init.b$Sigma <- matrix(rep(t(tmpSigma), init.b$K), ncol=dimb, byrow=TRUE)
            init.b$Li <- rep(init.b$Li, init.b$K)
          }else{
            if (length(init.b$Li) == prior.b$Kmax*LTb & prior.b$Kmax > init.b$K) init.b$Li <- init.b$Li[1:(init.b$K*LTb)]
            if (length(init.b$Li) != init.b$K*LTb) stop(paste("init", number, ".b$Li must be of length ", init.b$K*LTb, sep=""))
            init.b$Sigma <- matrix(NA, ncol=dimb, nrow=dimb*init.b$K)
            for (j in 1:init.b$K){
              tmpSigma <- matrix(0, nrow=dimb, ncol=dimb)
              tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- init.b$Li[((j-1)*LTb+1):(j*LTb)]
              tmpSigma <- tmpSigma %*% t(tmpSigma)
              err <- try(tmpSigma <- chol(tmpSigma), silent=TRUE)
              if (class(err) == "try-error") stop(paste("the ", j,"-th block of init", number, ".b$Li does not lead to a positive definite matrix", sep=""))
              tmpSigma <- chol2inv(tmpSigma)
              init.b$Sigma[((j-1)*dimb):(j*dimb),] <- tmpSigma
            }  
          } 
          rownames(init.b$Sigma) <- paste("j", rep(1:init.b$K, each=dimb), ".", rep(1:dimb, init.b$K), sep="")
          colnames(init.b$Sigma) <- paste("m", 1:dimb, sep="")        
          names(init.b$Li) <- paste("Li", rep(1:init.b$K, each=LTb), rep(naamLTb, init.b$K), sep="")
        }  
      }   
    }else{                   ### Sigma is checked and Li is computed from Sigma
      if (any(is.na(init.b$Sigma))) stop(paste("NA in init", number, ".b$Sigma"), sep="") 
      if (dimb == 1){
        if (length(init.b$Sigma) == 1) init.b$Sigma <- rep(init.b$Sigma, init.b$K)
        if (length(init.b$Sigma) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$Sigma <- init.b$Sigma[1:init.b$K]      
        if (length(init.b$Sigma) != init.b$K) stop(paste("init", number, ".b$Sigma must be of length ", init.b$K, sep=""))
        init.b$Sigma <- as.numeric(init.b$Sigma)
        names(init.b$Sigma) <- paste("Sigma", 1:init.b$K, sep="")      
        if (any(init.b$Sigma <= 0)) stop(paste("init", number, ".b$Sigma must be positive", sep=""))
        init.b$Li <- sqrt(1 / init.b$Sigma)
        names(init.b$Li) <- paste("Li", 1:init.b$K, sep="")      
      }else{
        if (!is.matrix(init.b$Sigma)) stop(paste("init", number, ".b$Sigma must be a matrix", sep=""))
        if (ncol(init.b$Sigma) != dimb) stop(paste("init", number, ".b$Sigma must have ", dimb, " columns", sep=""))
        if (nrow(init.b$Sigma) == dimb){
          if (any(init.b$Sigma[lower.tri(init.b$Sigma)] != t(init.b$Sigma)[lower.tri(init.b$Sigma)])) stop(paste("init", number, ".b$Sigma must be a symmetric matrix", sep=""))
          err <- try(Sigmainv <- chol(init.b$Sigma), silent=TRUE)
          if (class(err) == "try-error") stop(paste("Cholesky decomposition of init", number, ".b$Sigma failed", sep=""))
          Sigmainv <- chol2inv(Sigmainv)
          Litmp <- t(chol(Sigmainv))
          init.b$Li <- rep(Litmp[lower.tri(Litmp, diag=TRUE)], init.b$K)
        }else{
          if (nrow(init.b$Sigma) == prior.b$Kmax*dimb & prior.b$Kmax > init.b$K) init.b$Sigma <- init.b$Sigma[1:(init.b$K*dimb),]
          if (nrow(init.b$Sigma) != init.b$K*dimb) stop(paste("init", number, ".b$Sigma must have ", init.b$K, " times ", dimb, " rows", sep=""))
          init.b$Li <- numeric(0)
          for (j in 1:init.b$K){
            Sigmainv <- init.b$Sigma[((j-1)*dimb+1):(j*dimb),]
            if (any(Sigmainv[lower.tri(Sigmainv)] != t(Sigmainv)[lower.tri(Sigmainv)])) stop(paste(j, "-th block of init", number, ".b$Sigma is not symmetric", sep=""))
            err <- try(Sigmainv <- chol(Sigmainv), silent=TRUE)
            if (class(err) == "try-error") stop(paste("Cholesky decomposition of the ", j, "-th block of init", number, ".b$Sigma failed", sep=""))
            Sigmainv <- chol2inv(Sigmainv)
            Litmp <- t(chol(Sigmainv))
            init.b$Li <- c(init.b$Li, Litmp[lower.tri(Litmp, diag=TRUE)])
          }
        }
        rownames(init.b$Sigma) <- paste("j", rep(1:init.b$K, each=dimb), ".", rep(1:dimb, init.b$K), sep="")
        colnames(init.b$Sigma) <- paste("m", 1:dimb, sep="")              
        names(init.b$Li) <- paste("Li", rep(1:init.b$K, each=LTb), rep(naamLTb, init.b$K), sep="")
      }  
    }    
    
    ##### init.b:  gammaInv
    ##### ----------------------------------------------------  
    if (is.na(ib.gammaInv)){
      if (dimb == 1) init.b$gammaInv <- prior.b$zeta * bbVar
      else           init.b$gammaInv <- prior.b$zeta * diag(bbVar)
    }
    init.b$gammaInv <- as.numeric(init.b$gammaInv)
    if (length(init.b$gammaInv) == 1) init.b$gammaInv <- rep(init.b$gammaInv, dimb)
    if (length(init.b$gammaInv) != dimb) stop(paste("init", number, ".b$gammaInv must be of length ", dimb, sep=""))
    if (any(is.na(init.b$gammaInv))) stop(paste("NA in init", number, ".b$gammaInv"), sep="")
    if (any(init.b$gammaInv <= 0)) stop(paste("non-positive values in init", number, ".b$gammaInv"), sep="")    
    names(init.b$gammaInv) <- paste("gammaInv", 1:dimb, sep="")

    ##### init.b:  df
    ##### ----------------------------------------------------  
    if (is.na(ib.df)){
      if (prior.b$distribution == "MVT") init.b$df <- rgamma(init.b$K, shape=prior.b$gdf, rate=prior.b$hdf)
      else                               init.b$df <- rep(1000, init.b$K)
    }
    init.b$df <- as.numeric(init.b$df)
    if (length(init.b$df) == 1) init.b$df <- rep(init.b$df, init.b$K)
    if (length(init.b$df) != init.b$K) stop(paste("init", number, ".b$df must be of length ", init.b$K, sep=""))
    if (any(is.na(init.b$df))) stop(paste("NA in init", number, ".b$df"), sep="")
    if (any(init.b$df <= 0)) stop(paste("non-positive values in init", number, ".b$df"), sep="")        
    names(init.b$df) <- paste("df", 1:init.b$K, sep="")
    
    ##### init.b:  r
    ##### ----------------------------------------------------  
    if (is.na(ib.r)){
      if (dimb == 1){
        initz <- (init.b$b - scale.b$shift)/scale.b$scale
        MEANS <- matrix(rep(init.b$mu, I), ncol=init.b$K, byrow=TRUE)
        SDS   <- matrix(rep(sqrt(init.b$Sigma), I), ncol=init.b$K, byrow=TRUE)
        YY    <- matrix(rep(initz, init.b$K), ncol=init.b$K)
        WW    <- matrix(rep(init.b$w, I), ncol=init.b$K, byrow=TRUE)
        PROB  <- WW * dnorm(YY, mean=MEANS, sd=SDS)
      }else{
        initz <- (init.b$b - matrix(rep(scale.b$shift, I), ncol=dimb, byrow=TRUE))/matrix(rep(scale.b$scale, I), ncol=dimb, byrow=TRUE)
        PROB <- matrix(0, nrow=I, ncol=init.b$K)
        for (j in 1:init.b$K){
          MEANS <- init.b$mu[((j-1)*dimb+1):(j*dimb)]
          SIGMA <- init.b$Sigma[((j-1)*dimb+1):(j*dimb),]
          PROB[,j] <- init.b$w[j] * dMVN(initz, mean=MEANS, Sigma=SIGMA)        
        }        
      }
      sumPROB <- apply(PROB, 1, sum)
      sumPROB[sumPROB <= 0] <- 1
      PROB    <- PROB / matrix(rep(sumPROB, each=init.b$K), ncol=init.b$K, byrow=TRUE)
      init.b$r <- apply(PROB, 1, which.max)          
    }
    init.b$r <- as.numeric(init.b$r)
    if (length(init.b$r) != I) stop(paste("init", number, ".b$r must be of length ", I, sep=""))
    if (any(is.na(init.b$r))) stop(paste("NA in init", number, ".b$r", sep=""))
    if (any(init.b$r < 1) | any(init.b$r > init.b$K)) stop(paste("init", number, ".b$r out of the range (must lie between ", 1, " and ", init.b$K, ")", sep=""))
    names(init.b$r) <- paste("r", 1:I, sep="")    
  }else{
    init.b <- list(b=0, K=0, w=0, mu=0, Sigma=0, Li=0, gammaInv=0, df=0, r=0)
  }  
  
  return(init.b)
}

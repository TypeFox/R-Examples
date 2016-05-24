"ahaz.readold" <- function(surv, X, weights)
  {
    ## Purpose: (Old) function for reading survival data/covariates;
    ##          checks validity, format 
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Design matrix
    ##   weights    : Weight vector (nonnegative)
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    # Validity checks/coercions
    X<-as.matrix(X)
    if (!class(surv) == "Surv") stop("First argument must be a Surv object")
    if (!is.numeric(X)) stop("X must be a numeric matrix")
    if (nrow(surv) != nrow(X)) stop("Survival times/covariates have incorrect dimensions")
    if (nrow(surv) < 2) stop("too few observations")
    if (missing(weights)){
      weightsum <- nrow(X)
      weights <- rep(1, nrow(X)) / nrow(X)
    }
    else {
      if (length(weights) != nrow(surv))
        stop("Weights have incorrect dimension")
      if (sum(weights != 0) < 2 || any(weights < 0))
        stop("Negative weights or too few nonzero weights")
      weightsum <-sum(weights)
      weights   <- weights / weightsum
      
    }
    if (attr(surv,"type") != "right" && attr(surv,"type") != "counting")
      stop("Only right-censored or counting process data supported")
    if (any(duplicated(surv[surv[,ncol(surv)] == 1, ncol(surv)-1])))
        stop("Tied survival times not supported - break ties first")   
    if (any(is.na(X))) stop("Missing covariates not allowed") 
    if(attr(surv,"type") == "right"){right <- TRUE} else{right <- FALSE}

    # Number of observations/covariates
    nobs  <- as.integer(sum(weights>0))  
    nvars <- as.integer(ncol(X)) 

    names <- colnames(X)
    
    # If zero weights, remove observations
    if (any(weights == 0)) {
      X <- rbind(X[weights > 0,])
      surv <- surv[weights > 0,]
      weights <- weights[weights > 0]
    }

   # Left truncation requires two copies of X
    if (right) {
        ix       <- order(-surv[,1])  # Event times in decreasing order
        X        <- X[ix,]
        surv2    <- surv[ix,]
        times    <- surv2[,1]  # Exit times
        tdiff    <- c(-diff(surv2[,1]), min(surv2[,1]))  # Event time differences
        inout    <- weights[ix]   # Weighted at-risk-indicator
        wgt      <- abs(inout)
        death.yn <- surv2[,2]  # Event?
        atrisk   <- cumsum(inout)  # Weighted '#at risk'
        iatrisk  <- ifelse(atrisk==0,0,1/atrisk)
        tatrisk  <- times
      } else {
        eventtms <- c(surv[,1], surv[,2])
        ix <- order(-eventtms)
        death.yn <- c(rep(0, nobs), surv[,3])[ix]
        times    <- eventtms[ix]  # Entry/exit times
        tdiff    <- c(-diff(eventtms[ix]),0) # Event time differences
        inout    <- c(rep(c(-1, 1), each = nobs))[ix] * (rep(weights, 2)[ix])
        wgt      <- abs(inout/2)
        X        <- X[rep(1:nobs, 2)[ix],]
        atrisk   <- c(cumsum(inout)[-(2*nobs)],1) # Normalized 'number at risk'
        iatrisk  <- ifelse(atrisk==0,0,1/atrisk)
        tatrisk  <- rep(surv[,2]-surv[,1],2)[ix]
      }

    out <- list("X" = t(X), "times" = times,"tatrisk"=tatrisk, "tdiff" = tdiff,
                "death.yn" = death.yn, "ix" = ix, "surv" = surv,
                "inout" = inout, "wgt" = wgt,"iatrisk" = iatrisk, "nobs" = nobs,
                "nvars" = nvars, "weights" = weights*weightsum, "weightsum"=weightsum,
                "right" = right, "colnames" = names)    
    return(out)
  }


"ahaz.readnew" <- function(surv, X, weights, standardize = TRUE)
  {
    ## Purpose: Function for reading survival data/covariates;
    ##          checks validity, format 
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Design matrix
    ##   weights    : Weight vector (nonnegative)
    ##   standardize: Standardize X?
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    # Validity checks/coercions
    X<-as.matrix(X)
    if (!class(surv) == "Surv") stop("First argument must be a Surv object")
    if (!is.numeric(X)) stop("X must be a numeric matrix")
    if (nrow(surv) != nrow(X)) stop("Survival times/covariates have incorrect dimensions")
    if (nrow(surv) < 2) stop("too few observations")
    if (missing(weights)){
      weights <- rep(1, nrow(X))
    }
    else {
      if (length(weights) != nrow(surv))
        stop("Weights have incorrect dimension")
      if (sum(weights != 0) < 2 || any(weights < 0))
        stop("Negative weights or too few nonzero weights")
    }
    if (attr(surv,"type") != "right" && attr(surv,"type") != "counting")
      stop("Only right-censored or counting process data supported")
    if (any(duplicated(surv[surv[,ncol(surv)] == 1, ncol(surv)-1])))
        stop("Tied survival times not supported - break ties manually first")   
    if(attr(surv,"type") == "right"){right <- TRUE} else{right <- FALSE}

    # Number of observations/covariates
    nobs  <- as.integer(sum(weights>0))  
    nvars <- as.integer(ncol(X)) 

    # Variable names
    names <- colnames(X)
    
    # If any zero weights, remove these observations
    if (any(weights == 0)) {
      X <- rbind(X[weights > 0,])
      surv <- surv[weights > 0,]
      weights <- weights[weights > 0]
    }

    # Handle counting process data...
    
    # Copying X is an ugly and slow solution. Find a
    # suitable C-data structure instead.
    if(right){
      time1<-rep(0,nobs); time2<-surv[,1]; event<-surv[,2]
     } else {
      time1<-surv[,1]; time2<-surv[,2]; event<-surv[,3]; X<-rbind(X,X);
      weights<-rep(weights,2);
    }
    
    out<-list("X"=X,"time1"=time1,"time2"=time2,"event"=event,"surv"=surv,
                "weights"=weights,"standardize"=standardize,"colnames"=names,
              "nobs"=nobs,"nvars"=nvars,"rightcens"=right)
    return(out)
  }


"ahaz.predbackend"<-function(x, beta, type=c("residuals", "cumhaz"), colnames = NULL)
  {
    ## Purpose: Prediction backend used for both 'ahaz' and 'ahazpen' objects
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x       : list produced by ahaz.read
    ##   beta    : parameter estimate
    ##   type    : 'residuals' for ahaz residuals, 'cumhaz' for
    ##             Breslow cumulative hazard estimate
    ##   colnames: User-supplied column names (type='residuals' only)
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    type <- match.arg(type)
    
    if (missing(beta)) {stop("no beta specified")}
    if (!is.numeric(beta)) {stop("beta must be numeric")}
    if (length(beta) != nrow(x$X)) {stop("incorrect dimensions of X/beta")}
      
    a <- .C("ahbreslow",
            X       = as.double(x$X),
            tdiff   = as.numeric(x$tdiff),
            inout   = as.double(x$inout),
            iatrisk  = as.double(x$iatrisk),
            deathyn = as.integer(x$death.yn),
            n       = as.integer(length(x$tdiff)),
            p       = as.integer(x$nvars),
            beta    = as.numeric(beta),
            bres    = numeric(length(x$tdiff)),
            zbar    = numeric(length(x$tdiff) * x$nvars))
    
    bres <- a$bres
    zbar <- matrix(a$zbar, ncol=length(x$tdiff),nrow = x$nvars)
 
    if (type == "cumhaz") {
      if(x$right)
        out <- list("times" = rev(x$times),"event"=rev(x$death.yn), "cumhaz" = c(0,cumsum(rev(bres))),"br"=rev(bres))
      else
        out <- list("times" = rev(x$times)[-1],"event"=rev(x$death.yn)[-1], "cumhaz" = cumsum(rev(bres)))

      return(structure(out, class = "cumahaz"))

    } else
      {
        if (x$right)
          {
            start <- rep(0,x$nobs)
            times <- c(x$times,0)
          } else { 
            start <- x$surv[,1]
            times <- x$times
          }
        b <- .C("ahresid",
                start   = as.double(start),
                end     = as.double(x$surv[,ncol(x$surv)-1]),
                status  = as.double(x$surv[,ncol(x$surv)]),
                X       = as.double(x$X[,order(x$ix)[1:x$nobs]]),
                Zbar    = as.double(zbar),
                times   = as.double(times),
                tdiff   = as.double(x$tdiff),
                breslow = as.double(bres),
                beta    = as.double(beta),
                ntimes  = as.integer(length(x$tdiff)+1),
                p       = as.integer(x$nvars),
                nobs    = as.integer(x$nobs),
                resid   = as.double(matrix(0,nrow=x$nvars,ncol=x$nobs)),
                wgt     = as.double(x$weights))
        
        out<- t(matrix(b$resid, nrow = x$nvars,ncol = x$nobs))
        if(!is.null(colnames) && length(colnames) == x$nvars)
          colnames(out) <-colnames
        
        return(out)
      }
  }

"ahaz.linterp" <- function(lambda,s)
{
  ## Purpose: Interpolation for lambda sequences
  ## ---------------------------------------------------------------------
  ## This is the function 'lambda.interp' from  package
  ## 'glmnet' version 1.6, Jerome Friedman, Trevor Hastie, Rob Tibshirani.
  
  if(length(lambda)==1){# degenerate case of only one lambda
    nums=length(s)
    left=rep(1,nums)
    right=left
    sfrac=rep(1,nums)
  }
  else{
      s[s > max(lambda)] = max(lambda)
      s[s < min(lambda)] = min(lambda)
      k=length(lambda)
      sfrac <- (lambda[1]-s)/(lambda[1] - lambda[k])
      lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
      coord <- approx(lambda, seq(lambda), sfrac)$y
      left <- floor(coord)
      right <- ceiling(coord)
      sfrac=(sfrac-lambda[right])/(lambda[left] - lambda[right])
      sfrac[left==right]=1
    }
  return(list(left=left,right=right,frac=sfrac))
}


"ahaz.nzcoef" <- function(beta,bystep=FALSE){

  ## Purpose: Get nonzero part of a dgCMatrix
  ## ---------------------------------------------------------------------
  ## This is the function 'nzcoef' from package 'glmnet' version 1.6
  ## Jerome Friedman, Trevor Hastie, Rob Tibshirani.
  
  if(nrow(beta)==1){#degenerate case
    if(bystep)
      return(apply(beta,2,function(x)if(abs(x)>0)1 else NULL))
    else {if(any(abs(beta)>0))1 else NULL}
  }
  else{
    beta<-t(beta)
    which<-diff(beta@p)
    which<-seq(which)[which>0]
    if(bystep){
      nzel<-function(x,which)if(any(x))which[x]else NULL
      beta<-abs(as.matrix(beta[,which]))>0
      return(apply(beta,1,nzel,which))
    }
    else
      return(which)
  }
}

"ahaz.cvfolds" <- function (n, folds = 10)
  {
    ## Purpose: Return list of folds for cross validation
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   n    : number of observations
    ##   folds: number of folds
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    return(split(sample(1:n), rep(1:folds, length = n)))
  }

"error.bars" <- function(x, upper, lower, width = 0.01, ...)
{
    ## Purpose: Error bars for cross validation error curves
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x    : Vector of x-coordinates 
    ##   upper: Vector of upper limits
    ##   lower: Vector of lower limits
    ##   width: Width of error bars
    ##   ...  : Additional options for 'segments' function
    ## ----------------------------------------------------------------------
  
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
}

"ahaz.ginv"<-function(X)
  {
    ## Purpose: Generalized matrix inverse
    ## ----------------------------------------------------------------------
    ## Strongly inspired by the MASS implementation ginv
    svd<-svd(X)
    pos<-svd$d>sqrt(.Machine$double.eps)
    if(all(pos))
      return(svd$v%*%((1/svd$d)*t(svd$u)))
    else
      return(svd$v[,pos]%*%((1/svd$d[pos])*t(svd$u[,pos])))
  }

"ahaz.mse"<-function(x, beta)
  {
    ## Purpose: MSE from ahaz object and given beta
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'ahaz' object
    ##   beta  : beta estimate
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    if(length(beta)>1)
      return((t(beta) %*% x$D %*% beta - 2 * t(beta) %*% x$d))
    return(x$D*beta^2-2*beta*x$d)
  }

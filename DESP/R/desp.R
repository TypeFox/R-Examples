# DESP/R/desp.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (version 3) as published by
#  the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

desp <- 
function(X,lambda=0,gamma=0,settings=NULL){
  # estimation of the precision matrix Omega and of the matrix of outliers Theta 

  # read the sample size and the number of variables
  D <- dim(X);
  n <- D[1];               # n is the sample size
  p <- D[2];               # p is the dimension

  normE <- function(x){
    # Euclidean norm of a vector
    sqrt(sum(x^2))
  }

  norm1 <- function(x){
    # l1 norm of a vector
    sum(abs(x))
  }

  # settings
  lst = list("nThreads"=1, "thres"=1e-6, "maxIter"=2000, "diagElem"="RV", "OLS"=TRUE, "refit"=TRUE, "refit.min"=max(2, ceiling(n/2)), "OLS.min"=max(2, ceiling(n/2)), "sqrLasso"="CD", "sto"='0', "outParCorrB"=FALSE, "posDef"=TRUE, "SML.method"="MST", "SPT.spec"='1', "PML.thresh"=0.05, "PML.kappa"=0.05, "PML.tol"=1e-5, "rmSO"=TRUE, "rmSO.method"="MAD", "rmSO.mad.mult"=5, "rmSO.mad.constant"=1.4826, "rmSO.iqr.mult"=3)
  nSett <- names(lst)
  if(length(setdiff(names(settings),nSett))!=0) 
    {
    invalid <- paste0("'",setdiff(names(settings),nSett), sep="'", collapse = ' ')
    warning(gettextf("%s are not valid settings, hence ignored", invalid)) 
    }
  if (is.null(settings)) settings <- lst
  else {
    settings <- mapply(function(a,b) ifelse(is.null(b),a,b),lst, settings[nSett], SIMPLIFY=FALSE)
  }

  # available methods to solve the square-root Lasso
  mSqrL <- c(0,1)
  names(mSqrL) <- c("CD","SCS")

  # available methods to compute the diagonal entries of the precision matrix
  mDiag <- c("RV","AD","RML","SML","PML")

  # available methods to choose the path between two vertices
  mPath <- c("MST","SPT")
  # available variants for the SPT method
  vSPT <- c('1','2','2.1','3')

  # possible values for the parameter 'sto'
  vSto <- c('0','1','2')

  # available methods to detect simple outliers (with very large Euclidean norm)
  mSO <- c("MAD","IQR")

  # check the arguments of the function
  if(!is.matrix(X)){ stop("'X' must be a matrix");}
  if(!is.numeric(lambda) || lambda<0){ stop("'lambda' must be a non-negative scalar");}
  if(!is.numeric(gamma) || gamma<0){ stop("'gamma' must be a non-negative scalar");}
  if(!is.numeric(settings$nThreads) || settings$nThreads<1){ warning("'nThreads' must be an integer at least equal to 1, it has been set to default value 1"); settings$nThreads <- 1;}
  if(settings$nThreads!=1 && settings$sqrLasso == 'CD' && settings$sto=='1'){ warning("when using stochastic coordinate descent, 'nThreads' must be equal to 1, it has been set to default value 1"); settings$nThreads <- 1;}
  if(!is.numeric(settings$thres) || settings$thres<=0){ warning("'thres' must be a positive scalar, it has been set to default value 1e-6"); settings$thres <- 1e-6;}
  if(!is.numeric(settings$maxIter) || settings$maxIter<1){ warning("'maxIter' must be an integer at least equal to 1, it has been set to default value 2000"); settings$maxIter <- 2000;}
  if(!(settings$diagElem %in% mDiag)){ warning(gettextf("%s does not correspond to a valid method to compute the diagonal entries, %s is used instead", dQuote(settings$diagElem), dQuote('RV'))); settings$diagElem <- 'RV';}
  if(!is.logical(settings$OLS)){ warning("'OLS' must be a boolean, it has been set to default value TRUE"); settings$OLS <- TRUE;}
  if(!is.logical(settings$refit)){ warning("'refit' must be a boolean, it has been set to default value TRUE"); settings$refit <- TRUE;}
  if(settings$refit && (!is.numeric(settings$refit.min) || settings$refit.min<2)){ warning(gettextf("'refit.min' must be an integer at least equal to 2, it has been set to default value %d",max(2, ceiling(n/2)))); settings$refit.min <- max(2, ceiling(n/2));}
  if(settings$OLS && (!is.numeric(settings$OLS.min) || settings$OLS.min<2)){ warning(gettextf("'OLS.min' must be an integer at least equal to 2, it has been set to default value %d", max(2, ceiling(n/2)))); settings$OLS.min <- max(2, ceiling(n/2));}
  if(!(settings$sqrLasso %in% names(mSqrL))){ warning(gettextf("%s does not correspond to a valid method to solve the square-root Lasso, %s is used instead", dQuote(settings$sqrLasso), dQuote('CD'))); settings$sqrLasso <- 'CD';}
  if(settings$sqrLasso == 'CD' && !(settings$sto %in% vSto)){ warning(gettextf("%s does not correspond to a valid choice for the parameter 'sto', %s is used instead", dQuote(settings$sto), dQuote('0'))); settings$sto <- '0';}
  if(!is.logical(settings$outParCorrB)){ warning("'outParCorrB' must be a boolean, it has been set to default value FALSE"); settings$outParCorrB <- FALSE;}
  if(!is.logical(settings$posDef)){ warning("'posDef' must be a boolean, it has been set to default value TRUE"); settings$posDef <- TRUE;}
  if(settings$diagElem == 'SML' && !(settings$SML.method %in% mPath)){ warning(gettextf("%s does not correspond to a valid method to choose a path, %s is used instead", dQuote(settings$SML.method), dQuote('MST'))); settings$SML.method <- 'MST';}
  if(settings$diagElem == 'SML' && settings$SML.method == 'SPT' && !(settings$SPT.spec %in% vSPT)){ warning(gettextf("%s does not correspond to a valid choice to build the shortest path trees, %s is used instead", dQuote(settings$SPT.spec), dQuote('1'))); settings$SPT.spec <- '1';}
  if(settings$diagElem == 'PML' && (!is.numeric(settings$PML.thresh) || settings$PML.thresh<=0)){ warning("'PML.thresh' must be a positive scalar, it has been set to default value 0.05"); settings$PML.thresh <- 0.05;}
  if(settings$diagElem == 'PML' && (!is.numeric(settings$PML.kappa) || settings$PML.kappa<=0)){ warning("'PML.kappa' must be a positive scalar, it has been set to default value 0.05"); settings$PML.kappa <- 0.05;}
  if(settings$diagElem == 'PML' && (!is.numeric(settings$PML.tol) || settings$PML.tol<=0)){ warning("'PML.tol' must be a positive scalar, it has been set to default value 1e-5"); settings$PML.tol <- 1e-5;}
  if(!is.logical(settings$rmSO)){ warning("'rmSO' must be a boolean, it has been set to default value TRUE"); settings$rmSO <- TRUE;}
  if(settings$rmSO && !(settings$rmSO.method %in% mSO)){ warning(gettextf("%s does not correspond to a valid method to detect simple outliers, % is used instead", dQuote(settings$rmSO.method), dQuote('MAD'))); settings$rmSO.method <- 'MAD';}
  if(settings$rmSO && settings$rmSO.method == "MAD" && settings$rmSO.mad.mult<=0){ warning("'rmSO.mad.mult' must be a positive scalar, it has been set to default value 5"); settings$rmSO.mad.mult <- 5;}
  if(settings$rmSO && settings$rmSO.method == "MAD" && settings$rmSO.mad.constant<=0){ warning("'rmSO.mad.constant' must be a positive scalar, it has been set to default value 1.4826"); settings$rmSO.mad.constant <- 1.4826;}
  if(settings$rmSO && settings$rmSO.method == "IQR" && settings$rmSO.iqr.mult<=0){ warning("'rmSO.iqr.mult' must be a positive scalar, it has been set to default value 3"); settings$rmSO.iqr.mult <- 3;}

  if(settings$rmSO)
    {
    simpleOutliers <- DESP_SimpleOut(apply(X,1,normE), method=settings$rmSO.method, iqr.mult=settings$rmSO.iqr.mult, mad.mult=settings$rmSO.mad.mult, mad.constant=settings$rmSO.mad.constant)
    if(length(simpleOutliers)!=0)
      {
      X <- X[-simpleOutliers,]
      n <- n - length(simpleOutliers)
      }
    }

  # it is convenient to normalize X  
  X_n <- X/sqrt(n)
  u_n <- rep(1, n)/sqrt(n)
  # to avoid unnecessary recalculations
  Xc <- X_n-u_n%*%crossprod(u_n,X_n)

  # the default value of Theta is 0 (no outlier)
  Theta <- matrix(0,n,p)
  # some initialization for the matrix of the coefficients of regression B
  B <- diag(p)

  # the initial value of the objective function 
  obj <- sum(apply(Xc,2,normE)) # initializing B as the identity matrix
  
  # we determine the largest penalty parameter
  # above which all the entries of Theta will be zero
  lambda_max <- .Call("largestLambda",Xc)
  if (is.nan(lambda_max))
    {
    stop("NaN' occur when computing the largest value of 'lambda'")
    }

  if (lambda > lambda_max)
    {
    warning(gettextf("the penalty level %f is higher than %f, no outlier will be detected with such a high penalty", lambda, lambda_max))
    }

  if (lambda == 0 || lambda > lambda_max) 
    { # no outliers detection
    if (gamma == 0)
      { # estimate a dense precision matrix (in the moderate dimensional case for instance)
      B <- .Call("moderate_B", Xc, NULL)$B
      }
    else
      { # estimate a sparse precision matrix (in the high-dimensional case for instance)
      B <- .Call("srL_B", Xc, NULL, as.double(gamma),as.integer(mSqrL[settings$sqrLasso]),as.integer(settings$sto),as.integer(settings$nThreads))$B
      }
    }
  else
    {
	  # stopping criterion: initial value equal to the prescribed threshold
  	stop <- settings$thres
	  h <- 1
	  stopNA <- 0 # to remember when the loop exit on stop = NA (missing)
	  while(stop>=settings$thres && h <= settings$maxIter)
		  {
		  if (gamma == 0)
			  { # estimate a dense precision matrix (in the moderate dimensional case for instance)
			  Bt <- .Call("moderate_B", Xc, Theta - u_n %*% crossprod(u_n,Theta))
			  }
		  else
			  { # estimate a sparse precision matrix (in the high-dimensional case for instance) 
        Bt <- .Call("srL_B", Xc, Theta - u_n %*% crossprod(u_n,Theta), as.double(gamma),as.integer(mSqrL[settings$sqrLasso]),as.integer(settings$sto),as.integer(settings$nThreads)) 
			  }
			B <- Bt$B
			t <- Bt$t

		  Theta <- .Call("update_theta", Xc, B, t, lambda)

		  if (sum(is.nan(Theta))==0)
		    {
			  prec_obj <- obj
			  #obj <- sum(t)/2 + sum(t((crossprod(diag(n) - tcrossprod(u_n), X_n) %*% B - (diag(n) - tcrossprod(u_n)) %*% Theta)^2)/t)/2 + lambda*sum(apply(Theta,1,normE)) # faster, but use the transpose, which may cause a memory problem when the sample size is very large
			  obj <- sum(t)/2 + sum(apply((tcrossprod(tcrossprod(Xc,t(B))- (Theta - u_n %*% crossprod(u_n,Theta)), diag(1/sqrt(t))))^2, 1, sum))/2 + lambda*sum(apply(Theta,1,normE))
			  if (gamma != 0)
		            {
				  obj <- obj + gamma*sum(abs(B))
			    }
			  stop <- abs(prec_obj-obj)/obj
			  if (is.na(stop)){ stop <- -1; stopNA <- 1 }
		    }
      else{ stop <- -1 }
		  h <- h+1
		  }
		if (h>settings$maxIter){ warning(gettextf("'desp' has reached the maximum number of iterations, when it stopped, the ratio of the variation of the cost function was equal to %f", stop)) }
		if (stopNA){ warning("'desp' stopped as the stopping value is missing") }
    }


  if (sum(is.nan(Theta))!=0){ warning("NaN appear in computation, no outlier will be detected"); Theta <- matrix(0,n,p) }

  rowL1 <- apply(Theta,1,norm1)
  if (sum(rowL1>0) == n){ warning("all observations have been detected as outliers") }

  if(settings$refit && sum(rowL1>0) > 0 && sum(rowL1==0) >= settings$refit.min)
    {
    # remove the observations identified as outliers and reestimate the model
    settings$refit <- FALSE
    ref <- desp(X[rowL1==0,], lambda=0, gamma, settings)
    Omega <- ref$Omega
    mu <- ref$mu
    Psi <- ref$Psi
    B <- ref$B
    }
  else
    {
    SPC <- DESP_SqPartCorr(B,n)
    if(settings$OLS)
      {
      if (sum(rowL1==0) >= settings$OLS.min)
        {
        # the ordinary least squares are computed once the outliers removed
        B <- DESP_OLS_B(X[rowL1==0,] - rep(1,sum(rowL1==0)) %*% crossprod(rep(1,sum(rowL1==0)),X[rowL1==0,])/sum(rowL1==0),SPC)
        }
      }

    if(settings$diagElem == 'RV')
      {
      # estimation of the diagonal entries using residual variance
      phi_inv <- DESP_RV(X - u_n %*% crossprod(u_n,X),B,(Theta - u_n %*% crossprod(u_n,Theta))*sqrt(n))
      }
    else if(settings$diagElem == 'AD')
      {
      # estimation of the diagonal entries using average absolute deviation around the mean
      phi_inv <- DESP_AD(X - u_n %*% crossprod(u_n,X),B,(Theta - u_n %*% crossprod(u_n,Theta))*sqrt(n))
      }
    else if(settings$diagElem == 'RML')
      {
      # estimation of the diagonal of the precision matrix by likelihood maximization
      phi_inv <- DESP_RML(X - u_n %*% crossprod(u_n,X),B,(Theta - u_n %*% crossprod(u_n,Theta))*sqrt(n))
      }
    else if(settings$diagElem == 'SML')
      {
      if(settings$SML.method == 'MST')
        {
        # estimation of the diagonal of the precision matrix using minimum spanning trees
        phi_inv <- DESP_MST(X - u_n %*% crossprod(u_n,X),B,(Theta - u_n %*% crossprod(u_n,Theta))*sqrt(n))
        }
      else if(settings$SML.method == 'SPT')
        {
        # estimation of the diagonal of the precision matrix using shortest path trees
        phi_inv <- DESP_SPT(X - u_n %*% crossprod(u_n,X),B,settings$SPT.spec,(Theta - u_n %*% crossprod(u_n,Theta))*sqrt(n))
        }
      }
    else if(settings$diagElem == 'PML')
      {
      # estimation of the diagonal of the precision matrix by penalized likelihood minimization
      phi_inv <- DESP_PML(X - u_n %*% crossprod(u_n,X),B,settings$PML.thresh,settings$PML.kappa,settings$PML.tol,(Theta - u_n %*% crossprod(u_n,Theta))*sqrt(n))
      }

    Psi <- diag(p)
    Psi[SPC>0] <- sqrt(SPC[SPC>0])*sign(B[SPC>0])
    Omega <- Psi * tcrossprod(sqrt(phi_inv))

    if(settings$posDef)
      {
    	# test for positive definiteness
    	cholStatus <- try(chol(Omega), silent = TRUE)
    	cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    	if (cholError)
        {
	      Omega <- Matrix::nearPD(Omega)$mat
        if(settings$outParCorrB)
          {
	        Psi <- Matrix::nearPD(Psi)$mat
          B <- tcrossprod(Omega, diag(1/phi_inv))
          }
    	  }
      }
    mu <- crossprod(X_n - Theta %*% MASS::ginv(B), u_n)

    }

  if(settings$outParCorrB)
    {
    res <- list(Omega=Omega,Theta=Theta,mu=mu,Psi=Psi,B=B)
    }
  else
    {
    res <- list(Omega=Omega,Theta=Theta,mu=mu)
    }
  if(settings$rmSO)
    {
    res$simpleOut <- simpleOutliers
    }
  class(res) <- "desp"
  return(res)
}




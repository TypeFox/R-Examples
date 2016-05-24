predict.monfd <- function(object, newdata=NULL, Lfdobj=0,
                            returnMatrix=FALSE, ...) {

  if(is.null(newdata))newdata <- object$argvals
##
## 1.  eval.monfd
##
  evalMon <- eval.monfd(newdata, object$Wfdobj, Lfdobj, returnMatrix)
##
## 2.  beta
##
  beta <- object$beta
  {
    if(length(dim(beta))<2){
      if(length(dim(evalMon))<2){
        be <- beta[2]*evalMon
        if(Lfdobj<1)
          be <- beta[1]+be
        return(be)
      }
      else
        stop('beta does not match eval.monfd(...)')
    }
    else {
      nem <- dim(evalMon)
      if(length(dim(beta)<3)) {
        if(length(nem)==2){
          be <- (evalMon*rep(beta[2,], each=nem[1]))
          if(Lfdobj<1)
            be <- (be+rep(beta[1,], each=nem[1]))
          return(be)
        }
        else
          stop('beta does not match eval.monfd(...)')
      }
      else {
        if(length(nem)==3){
          be <- (evalMon*rep(beta[2,,], each=nem[1]))
          if(Lfdobj<1)
            be <- (be+rep(beta[1,,], each=nem[1]))
          return(be)
        }
        else
          stop('beta does not match eval.monfd(...)')
      }
    }
  }

}

fitted.monfd <- function(object, ...){
  predict(object)
}

residuals.monfd <- function(object, ...){
  pred <- predict(object)
  object$y-pred
}

eval.monfd <- function(evalarg, Wfdobj, Lfdobj=int2Lfd(0), returnMatrix=FALSE) {
  #  Evaluates a monotone functional data observation, or the value of a linear
  #  differential operator LFD applied to the object,
  #  at argument values in an array EVALARGS.
  #  Functional data object LFD, if an integer, defines NDERIV, the
  #  order of derivative to be evaluated.
  #  Functional data object LFD, if a fd object, defines weight
  #  functions for computing the value of a linear differential operator
  #  applied to the functions that are evaluated.

  #  A monotone functional data object h  is in the form

  #           h(x) = [D^{-1} exp Wfdobj](x)

  #  where  D^{-1} means taking the indefinite integral.
  #  The interval over which the integration takes places is defined in
  #  the basisfd object in WFD.

  #  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
  #               from a call to function BsplineS.  See this function for
  #               enabling this option.

  #  Last modified 9 May 2012  by Jim Ramsay
  #  check Wfdobj

  if (!inherits(Wfdobj, "fd")) stop("Wfdobj is not a fd object.")

  #  extract number of variables and curves from coefficient matrix for Wfdobj

  coef  <- Wfdobj$coefs
  if (is.vector(coef)) coef <- as.matrix(coef)
  coefd <- dim(coef)
  ndim  <- length(coefd)
  if (ndim == 2) {
    ncurve <- coefd[2]
    nvar   <- 1
  } else {
    ncurve <- coefd[2]
    nvar   <- coefd[3]
  }

  #  determine if LFDOBJ is an integer

  Lfdobj <- int2Lfd(Lfdobj)

  if (!is.integerLfd(Lfdobj)) stop(
		"LFDOBJ is not an integer operator.")

  nderiv <- Lfdobj$nderiv

  n  <- length(evalarg)

  hmat <- array(0,c(n,ncurve,nvar))

  if (nderiv >= 2) Dwmat  <- getbasismatrix(evalarg, Wfdobj$basis, 1,
                                            returnMatrix)
  if (nderiv == 3) D2wmat <- getbasismatrix(evalarg, Wfdobj$basis, 2,
                                            returnMatrix)
  basislist = vector("list", 15)
  
  for (ivar in 1:nvar) {
    for (icurve in 1:ncurve) {

  	if (nderiv == 0) {
    	  if (ndim == 2) {
          if (ncurve == 1) {
            hmat[,icurve,ivar] <- monfn(evalarg, Wfdobj, basislist, 
                                        returnMatrix)
          } else {
            hmat[,icurve,ivar] <- monfn(evalarg, Wfdobj[icurve], basislist, 
                                        returnMatrix)
          }
        } else {
            hmat[,icurve,ivar] <- monfn(evalarg, Wfdobj[icurve,ivar], basislist,
                                        returnMatrix)
        }
  	}

  	if (nderiv == 1) {
    	  if (ndim == 2) {
            hmat[,icurve,ivar] <- exp(eval.fd(evalarg, Wfdobj[icurve], 0,
                                            returnMatrix))
        } else {
            hmat[,icurve,ivar] <- exp(eval.fd(evalarg, Wfdobj[icurve,ivar], 0,
                                            returnMatrix))
        }
  	}

  	if (nderiv == 2) {
        if (ndim == 2) {
          temp = (Dwmat %*% coef[,icurve])*
                                 exp(eval.fd(evalarg, Wfdobj[icurve],  0,
                                            returnMatrix))
    	    hmat[,icurve,ivar] <- as.vector(temp)
        } else {
          temp = (Dwmat %*% coef[,icurve])*
                                 exp(eval.fd(evalarg, Wfdobj[icurve,ivar], 0,
                                            returnMatrix))
    	    hmat[,icurve,ivar] <- as.vector(temp)
        }
  	}

  	if (nderiv == 3) {
        if (ndim == 2) {
    	    hmat[,icurve,ivar] <- as.vector(((D2wmat %*% coef[,icurve]) +
                                 (Dwmat  %*% coef[,icurve])^2)*
                                  exp(eval.fd(evalarg, Wfdobj[icurve], 0,
                                            returnMatrix)))
        } else {
    	    hmat[,icurve,ivar] <- as.vector(((D2wmat %*% coef[,icurve,ivar]) +
                                 (Dwmat  %*% coef[,icurve,ivar])^2)*
                                  exp(eval.fd(evalarg, Wfdobj[icurve,ivar],
                                            0, returnMatrix)))
        }
  	}

  	if (nderiv > 3) stop ("Derivatives higher than 3 are not implemented.")

    }
  }

  if (nvar == 1) hmat <- as.matrix(hmat[,,1])

  if((!returnMatrix) && (length(dim(hmat)) == 2)){
    return(as.matrix(hmat))
  }
  return(hmat)

}

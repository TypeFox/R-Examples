#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## regress.lars:
##
## fit y2 ~ y1  using lasso.


'regress.lars' <-
  function(y1, y2, method="lasso", validation="CV", verb=0)
{
  ## number of regressions
  numreg <- ncol(y2); if(is.null(numreg)) numreg <- 1
    
  ## number of observatios and predictors
  numobs <- nrow(y1); if(is.null(numobs)) numobs <- length(y1)
  numpred <- ncol(y1); if(is.null(numpred)) numpred <- 1

  ## decide on cross validation method depending on number of data points
  if(validation == "Cp" && numobs <= numpred) validation <- "CV"
  if(numobs <= 10) validation <- "LOO";
  if(validation == "LOO") K <- numobs
  else K <- 10
  actual.method <- paste(method, "-", validation, sep="")

  ## add to progress meter
  if(verb > 0) cat(paste("using ", method, " (", validation, ") ", sep=""))

  ## choose lambda (could be different for each response)
  nzb <- rep(NA, numreg)
  if(method != "lars" && method != "stepwise") lambda <- rep(NA, numreg)
  else lambda <- NULL
  bvec <- matrix(NA, nrow=numpred+1, ncol=numreg)
  res <- matrix(NA, nrow=numobs, ncol=numreg)
  if(verb > 0) cat("ncomp:")

  ## don't use gram matrix when m > 500
  use.Gram <- TRUE
  if(numpred > 500) use.Gram <- FALSE

  if(numreg == 1) y2 <- matrix(y2, ncol=numreg)
  for(i in 1:numreg) {

    if(validation == "Cp") {

      ## use (Mallows) Cp method
      reglst <- lars(x=y1, y=y2[,i],type=method,intercept=TRUE, use.Gram=use.Gram)
      if(sum(is.nan(reglst$Cp))) ## unsuccessful
        return(y1, y2, method=method, validation="CV", verb=verb)
      s <- as.numeric(names(which.min(reglst$Cp))) + 1
      mode <- "step"

    } else {
      
      ## use Cross-Validation (possibly LOO)      
      ## num.fractions could be passed in same as ncomp.max
      cv <- cv.lars(x=y1,y=y2[,i],type=method,K=K,intercept=TRUE,use.Gram=use.Gram,
                    plot.it=FALSE)
      
      ## choose with with "one-standard error rule"
      wm <- which.min(cv$cv)
      tf <- cv$cv < cv$cv[wm] + cv$cv.error[wm]
      s <- cv$index[(1:100)[tf][1]]
      ## s <- cv$fraction[(1:100)[tf][1]]
      mode <- "fraction"
    
      ## get the lasso fit with fraction f
      reglst <- lars(x=y1,y=y2[,i],type=method,intercept=TRUE,use.Gram=use.Gram)
      
    }

    ## extract the coefficients
    co <- coef(reglst, s=s, mode=mode)
  
    ## extract regression coefficients and parameters
    y1co <- drop(y1 %*% co)
    icept <- reglst$mu - mean(y1co)
    bvec[,i] <- c(icept, co)
    
    ## use predict to get the residuals
    res[,i] <- y2[,i] - (icept + y1co)

    ## calculate the number of nonzero coefficients, and extract lambda
    nzb[i] <- sum(co != 0)
    if(!is.null(lambda)) lambda[i] <- get.lambda(reglst, s, mode)
    if(verb > 0) cat(paste(" ", nzb[i], sep=""))

  }

  ## possibly cap off the print the number of components used
  if(verb > 0) cat(paste(" of ", min(numobs, numpred), sep=""))

  return(list(method=actual.method, ncomp=nzb, b=bvec, res=res, lambda=lambda))
}


## get.lambda:
##
## function to extract the lambda penalization parameter
## from the lars object (obj) corresponding to the
## fraction setting (s).  No checking of inputs is
## currently provided

get.lambda <- function(obj, s, mode)
  {
    ## base case
    if(s == 0) return(max(obj$lambda))
    
    ## extract the regression coefficients for the desired fraction
    beta <- coef(obj, s=s, mode=mode)

    ## get the non-zero entries
    bnz <- beta != 0
    b <- beta[bnz]

    ## extract the norm of OLS coefficients
    bls2 <- sum(coef(obj, s=0, mode="lambda")^2)#[nz]

    ## objective function to be minimized for x
    of <- function(x, obj, bnz, b, bls2)
      {
        beta <- coef(obj, s=x, mode="lambda")
        zdiff2 <- sum((bnz - (beta != 0))^2)
        bdiff2 <- sum((b - beta[bnz])^2)/bls2
        r <- bdiff2 + zdiff2
        return(r)
      }

    ## range of lambda values to optimize over
    lrange <- c(0, max(obj$lambda))

    ## find the corresponding lambda
    r <- optimize(of, interval=lrange, obj=obj, bnz=bnz,
                  b=b, bls=bls2)$minimum

    ## return the lambda that was found
    return(r)
  }

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


## opt.ridge:
##
## function which determines what to optimize in order to
## chose lambda -- LOO or GCV can be used.  Generally LOO
## is used for p >= n regressions

opt.ridge <- function(lambda, y2, y1, LOO=FALSE)
  {
    ## sanity check
    if(lambda <= 0) return(Inf)

    ## using Leave-One-Out Cross Validation
    if(LOO) {
      res <- rep(NA, length(y2))
      for (omit in 1:length(y2)) {
        x <- y1[-omit, ,drop=FALSE]; y <- y2[-omit]
        beta <- coef(lm.ridge(y~x, lambda=lambda))
        fit <- y1[omit, ,drop = FALSE] %*% (beta[-1]) + beta[1]
        res[omit] <- (y2[omit] - fit)^2
      }
      r <- mean(res)
      r <- r + 2*sd(res)
    }

    ## using GCV
    else r <- lm.ridge(y2 ~ y1, lambda=lambda)$GCV

    ## print(c(ncol(y1), length(y2), lambda, r))
    return(r)
  }



## regress.ridge:
##
## fit y2 ~ y1  using ridge regression
##
## NOTE: this does not seem to work for the regression with
## "big-p small-n" -- it essentially gives zero-residulas
## regardless of the method used to choose the ridge
## constant (lambda)

'regress.ridge' <-
  function(y1, y2, verb=0)
{
  ## number of regressions
  numreg <- ncol(y2); if(is.null(numreg)) numreg <- 1
    
  ## number of observatios and predictors
  numobs <- nrow(y1); if(is.null(numobs)) numobs <- length(y1)
  numpred <- ncol(y1); if(is.null(numpred)) numpred <- 1

  ## add to progress meter
  method <- "ridge"
  if(verb > 0) cat(paste("using ", method, " (LOO) ", sep=""))

  ## choose lambda (could be different for each response)
  lambda <- rep(NA, numreg)
  bvec <- matrix(NA, nrow=numpred+1, ncol=numreg)
  res <- matrix(NA, nrow=numobs, ncol=numreg)
  if(verb > 0) cat("lambda:")

  ## a hueristic for lower-bounding and upper-bounding of lambda
  lmin <- 0
  if(numpred > numobs) lmin <- log(numpred-numobs)
  lower <- lmin; upper <- numpred
  
  if(numreg == 1) y2 <- matrix(y2, ncol=numreg)
  for(i in 1:numreg) {

    ## actually get the ridge reggression constant (lambda)
    while(1) { ## loop out bounds of the optimization window
      lam <- optimize(opt.ridge, lower=lower, upper=upper,
                      y2=y2[,i], y1=y1)$minimum
      if(round(lam) != upper) break;
      lower <- lam; upper <- 2*upper
    }

    ## use that lambda, and re-fit the model to get the coefficients
    reglst <- lm.ridge(y2[,i] ~ y1, lambda=lam)

    ## get the coefficients and residulas and check that they
    ## have led to a valid regression
    ## for(j in 1:2) {

      ## get the regression coefficients and lambda estimate
      bvec[,i] <- coef(reglst)
      
      ## check that the residuals aren't too small
      res[,i] <- y2[,i] - (bvec[1,i] + y1 %*% bvec[-1,i])
      ##if(any(res[,i]^2 < .Machine$double.eps))
        ##{ reglst <- lm.ridge(y2[,i] ~ y1, lambda=Inf); next; }
      
      ## get lambda
      lambda[i] <- reglst$lambda
      if(verb > 0) cat(paste(" ", signif(lambda[i],5), sep=""))
      ##break;
    ##}
  }

  ## print lambda range
  if(verb > 0) cat(paste(" in range [", signif(lmin,5),
                         ",", upper, "]", sep=""))

  return(list(method=method, ncomp=lambda, lambda=lambda, b=bvec, res=res))
}

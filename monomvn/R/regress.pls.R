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


## regress.pls:
##
## fit y2 ~ y1  using pcr or plsr

`regress.pls` <-
function(y1, y2, method="plsr", ncomp.max=Inf, validation="CV", 
         verb=0, quiet=TRUE)
  {
    ## number of regressions
    numreg <- ncol(y2); if(is.null(numreg)) numreg <- 1
    
    ## number of observatios and predictors
    numobs <- nrow(y1); if(is.null(numobs)) numobs <- length(y1)
    numpred <- ncol(y1); if(is.null(numpred)) numpred <- 1

    ## force plsr when there are rediculiously few data points
    if(numobs <= 3 && method == "pcr") {
      if(!quiet) warning("forcing plsr when fewer than 4 non-missing entries")
      method <- "plsr"
    }
    
    ## possibly change cross validation method depending on number of data points
    if(numobs <= 10) validation <- "LOO"
    actual.method <- paste(method, "-", validation, sep="")
    
    ## add to progress meter
    if(verb > 0) cat(paste("using ", method, " (", validation, ") ", sep=""))
    
    ## maximum number of components to include in the regression
    ncomp.max <- min(ncomp.max, numpred, numobs-1)
    
    ## use principal least squares or principal component regression
    suppressWarnings({
      if(method == "plsr") reglst <- plsr(y2~y1,validation=validation,ncomp=ncomp.max)
      else if(method == "pcr") reglst <- pcr(y2~y1,validation=validation,ncomp=ncomp.max)
      else stop(paste("regression method", method, "unknown"))
    })
    
    ## choose number of components (could be different for each response)
    ncomp <- rep(NA, numreg)
    bvec <- matrix(NA, nrow=numpred+1, ncol=numreg)
    res <- matrix(NA, nrow=numobs, ncol=numreg)
    a <- matrix(RMSEP(reglst, estimate="CV", intercept=FALSE)$val, nrow=numreg)
    a[!is.finite(a)] <- 1
    if(verb > 0) cat("ncomp:")
    for(i in 1:numreg) {
      ## in future, the "one standard error" rule should be used here instead
      ## but the se=TRUE argument was not implemented yet in RMSEP
      ncomp[i] <- which.min(a[i,] + (max(a[i,])-min(a[i,]))*seq(0,1,length=ncol(a)))
      ## plot(a[i,]); abline(v=ncomp[i]); readline("pr: ")
      if(verb > 0) cat(paste(" ", ncomp[i], sep=""))
      bvec[,i] <- matrix(coef(reglst, ncomp=ncomp[i], intercept=TRUE), ncol=numreg)[,i]
      res[,i] <- reglst$resid[,i,ncomp[i]]
    }
    if(verb > 0) cat(paste(" of ", ncomp.max, sep=""))
    
    ## return method, mean vector, and mean-squared of residuals
    return(list(method=actual.method, ncomp=ncomp, bvec=bvec, res=res))
  }


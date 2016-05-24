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


## regress:
##
## fit y ~ x  using a linear model.
##
## This is basically the switch function for the regressions
## used by the monomvn algorithm.
## If p*nrow(x) >= p*ncol(x) then use principal least squares
## (pls) or principal component (pc) regression, ridge regression
## or one of the lars methods -- instead of lsfit

`regress` <-
function(X, y, method=c("lsr", "plsr", "pcr", "lasso", "lar",
                 "forward.stagewise", "stepwise", "ridge", "factor"),
         p=0.0, ncomp.max=Inf, validation=c("CV", "LOO", "Cp"),
         verb=0, quiet=TRUE)
  {
    ## save the call
    cl <- match.call()
    
    ## check method argument
    method <- match.arg(method)
    if(method == "lsr") p <- 1
    
    ## check p argument, whose usage varies depending on method="factor"
    if(method != "factor") {
      if(length(p) != 1 || p > 1 || p < 0) {
        warning("should have scalar 0 <= p <= 1, using p=1")
        p <- 1
      }
    } else { ## method == "factor"
      if(length(p) != 1 || p < 1) {
        warning("when method=\"factor\" pust have p >= 1, using  p=1")
        p <- 1
      }

      ## save p in num.facts and set p=0
      num.facts <- min(p, ncol(X));  p <- 0
    }

    ## check ncomp.max argument
    if(length(ncomp.max) != 1 || ncomp.max < 1) {
      warning("should have integer 1 <= ncomp.max, using default ncomp=Inf")
      ncomp.max <- Inf
    }
    
    ## check validation argument
    validation <- match.arg(validation)
    if(validation == "Cp" && (method == "plsr" || method == "pcr"))
      stop("Cp model selection is only valid for lars models, not plsr or pcr")

    ## check verb argument
    if(length(verb) != 1 || verb < 0)
      stop("verb should be a positive scalar")

    ## check quiet argument
    if(length(quiet) != 1 || !is.logical(quiet))
      stop("quiet should be a logical scalar")

    ## done checking the arguments
    
    ## number of observatios and predictors
    numobs <- nrow(X); if(is.null(numobs)) numobs <- length(X)
    numpred <- ncol(X); if(is.null(numpred)) numpred <- 1

    ## start progress meter
    if(verb > 0) cat(paste(numobs, "x", numpred, " ", sep=""))

    ## use non-LS regression when usepler-times the number of columns
    ## in the regression is >= the number of rows (length of y)
    if(!is.null(dim(X)) && (numpred > 1) && (numpred >= p*numobs)) {

      ## choose the regression method
      if(method == "plsr" || method == "pcr")
        ret <- regress.pls(X, y, method, ncomp.max, validation, verb, quiet)
      else if(method == "ridge") ret <- regress.ridge(X, y, verb)
      else if(method == "factor") ret <- regress.fact(X, y, num.facts, verb)
      else ret <- regress.lars(X, y, method, validation, verb)

      ## least squares regression
    } else ret <- regress.ls(X, y, verb)

    ## calculate the mean-square residuals
    ## S <- (numobs-1)*cov(ret$res)/numobs
    S <- cov(ret$res)

    ## make sure the regression was non-signuar.  If so,
    ## use force a pls (or maybe lars or ridge) regression & print a warning
    if(sum(S) == 0) {
    ##mres <- apply(S, 2, function(x){ mean(abs(x)) })
    ##if(any(mres < 1e-7)) {

      ## nothing to do if lar was not used
      if(p == 0 || ret$method != "lsr")
        stop(paste("singular", method, "regression"))
      
      if(!quiet) ## warn if lsr
        warning(paste("singular least-squares ", nrow(X), "x", numpred,
                    " regression, forcing ", method, " ", sep=""))
      if(verb > 0) cat("[FAILED], ")

      ## force a parsimonious regression
      return(regress(X, y, method, p=0, verb=verb))
    }

    ## sanity check -- could hallen if all X's are the same
    if(any(!is.finite(ret$b)) || any(!is.finite(S)))
      stop("encountered NA in regression")

    ## return method, mean vector, and mean-squared of residuals
    return(list(call=cl, method=ret$method, ncomp=ret$ncomp,
                lambda=ret$lambda, b=ret$b, S=S))
  }

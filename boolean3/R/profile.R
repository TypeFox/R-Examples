### ----------------------------------------------------------------------------
### This file is part of boolean3
###
### Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
###
### boolean3 represents a substantial re-write of the original boolean package
### developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
### was developed under the direction of Bear Braumoeller and with support from
### The Ohio State University's College of Social and Behavioral Sciences.
###
### boolean3 is free software: you can redistribute it and/or modify it under
### the terms of the GNU General Public License as published by the Free
### Software Foundation, either version 3 of the License, or (at your option)
### any later version.
###
### This program is distributed in the hope that it will be useful, but WITHOUT
### ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
### FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
### more details.
###
### You should have received a copy of the GNU General Public License along with
### this program.  If not, see <http://www.gnu.org/licenses/>.
###
### ----------------------------------------------------------------------------

##' Calculate estimated likelihood-profiles.
##'
##' This function calculates log-likelihood profiles for the selected
##' variables. Despite the function name, these are not true profile likelihoods
##' as they hold all other coefficients fixed at their MLE.
##'
##' @param obj object of \code{boolean-class} containing a fit boolean model.
##' @param method estimation method to use
##' @param vars numeric vector selecting a set of covariates from the fitted
##' model  
##' @param k integer indicating the number of points at which the
##' log-likelihood should be calculated. 
##' @param as.table logical (default \code{TRUE}), to be passed to
##' \code{xyplot}.
##' @param scales list of settings for the scales argument passed to
##' \code{xyplot}. 
##' @param between numeric specifying the space between panels.
##' @param main string, plot title
##' @param xlab string, the \code{x}-axis label. 
##' @param ylab string, the \code{y}-axis label. 
##' @param ... Additional arguments to pass to \code{\link{xyplot}}. See that
##' cumentation for details.
##' @return Returns an object of \code{\link{boolprof-class}}, the default
##' action being to present the default plot.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @examples
##' \dontrun{
##'
##' ## Note: This example assumes a boolean model has already been fit.
##'
##' ## Display the contours of the likelihood given a change the value of
##' ## the coefficients. 
##' (prof <- boolprof(fit))
##' 
##' ## Extract the plots for x1_a and x4_b.
##' plot(prof, y = c("x1_a", "x4_b"))
##' plot(prof, y = c(1, 3), scales = list(y = list(relation = "free")))
##' 
##' ## You can also use variable or index matching with boolprof to select 
##' ## particular covariates of interest.
##' boolprof(fit, vars = c(1, 3))
##' boolprof(fit, vars = c("x1_a", "x4_b"))
##' }
##' @export
boolprof <- function(obj, method=names(obj$model.fit)[1], vars=1:obj$k, k=50,
                     as.table=TRUE, scales=list(x=list(relation="free")),
                     between=list(x=1, y=1), main="Estimated likelihood profiles",
                     xlab="beta", ylab="Log-likelihood", ...) {

  ## get index values of vars if names are supplied.
  if (is.character(vars))
    vars <- match(vars, obj$coef.labels)
  
  ## Get summary object and number of covariates to be calculated.
  smry <- summary(obj)$methods[[method]]
  v <- length(vars)

  ## Set up the default profile matrix, which is a matrix of k rows containing
  ## the estimated coefficient values.
  prof <- matrix(smry$coef[,1], nrow=k, ncol=length(smry$coef[,1]),
                 byrow=TRUE)

  ## Set up the results storage and index matrices.
  results <- matrix(0, nrow=k*length(vars), ncol=3)
  idx  <- cbind(seq(1, k*v, by=k), seq(k, k*v, by=k))
  
  ## Iterate through the variables, calculating the profile of the
  ## log-likelihood centered at the value of the covariate of interest.
  for (i in 1:v) {
    j <- vars[i]
      
    ## Get est coefficient and se for the variable of interest.
    b  <- smry$coef[j, 1]
    se <- smry$coef[j, 2]
    
    ## Set points at which to calculate the log-likelihood; slightly wider than
    ## 2x the se on each side.
    x <- seq(b-(se*2.05), b+(se*2.05), length=k)
    prof[,j] <- x
    
    ## Calculate loglik at all points. Remember, calc_lik returns the negative
    ## of the log-likelihood.
    llik <- -1 * apply(prof, 1, obj$calc_lik)

    ## Return the prof matrix to its original state.
    prof[,j] <- b

    ## Create the results matrix and place it in the storage matrix.
    results[idx[i,1]:idx[i,2],] <- cbind(llik, x, j)
  }

  results <- data.frame(results)
  results[,3] <- factor(obj$coef.labels[results[,3]], levels=obj$coef.labels)
  names(results) <- c("llik", "x", "var")
  
  plot <- xyplot(llik ~ x | var, data=results, type="l", as.table=as.table,
                 scales=scales, between=between, main=main, xlab=xlab,
                 ylab=ylab, ...)

  ## Construct boolprof object and return.
  profobj <- list(est=results, coef.labels=obj$coef.labels,
                  default.plot=plot)
  class(profobj) <- c("boolprof", class(profobj))
  return(profobj)
}

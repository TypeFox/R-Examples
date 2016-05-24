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


##' Calculate predicted probabilities
##'
##' This function calculates predicted probabilities for the selected covariate
##' profiles.
##'
##' @param obj object of \code{boolean-class} containing a fit boolean model.
##' @param vars vector selecting a set of covariates from the fitted
##' model. This can be a character vector of covariate names (as output from
##' \code{summary(obj)}), or a numeric vector indexing the desired covariates.
##' @param newdata data.frame with the same structure as model.matrix(boolean). 
##' @param k integer indicating the number of points at which the
##' predicted probability should be calculated. 
##' @param conf.int logical; should confidence intervals be simulated.
##' @param n number of draws to take from the estimated parameter space.
##' @param as.table logical (default \code{TRUE}), to be passed to
##' \code{xyplot}.
##' @param scales list of settings for the scales argument passed to
##' \code{xyplot}. 
##' @param between numeric specifying the space between panels.
##' @param xlab string, the \code{x}-axis label. 
##' @param ylab string, the \code{y}-axis label. 
##' @param ... Additional arguments to pass to \code{\link{xyplot}}. See that
##' documentation for details.
##' @return Returns an object of \code{\link{boolprob-class}}, the default
##' action being to present the default plot.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @examples
##' \dontrun{
##'
##' ## Note: This example assumes a boolean model has already been fit.
##'
##' ## Plot predicted probabilities for a fitted model. Either vars or 
##' ## newdata *must* be specified.
##' boolprob(fit, vars = c("x1_a", "x4_b"))
##' boolprob(fit, vars = c(2, 3, 4, 6))
##' 
##' ## Specifying conf.int = TRUE produces simulated confidence intervals. 
##' ## The number of samples to pull from the distribution of the estimated
##' ## coefficients is controlled by n; n=100 is default. This can take a
##' ## while.
##' (prob <- boolprob(fit, vars = c(2, 3, 4, 6), n = 1000, conf.int = TRUE))
##' }
##' @export
boolprob <- function(obj, vars=NULL, newdata=NULL, k=50, conf.int=FALSE, n=100,
                     as.table=TRUE, scales=list(x=list(relation="free")),
                     between=list(x=1, y=1), xlab="x",
                     ylab="Predicted probability", ...)
{
    
    ## Define data points at which predicted probabilities should be
    ## calculated. If newdata is supplied, it is retained as-is (though placed
    ## into a list).
    newdata <- set_newdata(obj, vars, newdata, k)
    
    ## Calculate the predicted probabilities for the covariate profiles defined
    ## in newdata.
    prob <- calc_prob(newdata, conf.int, n, obj)
    
    ## Construct plot.
    if (conf.int) {
        plot <- xyplot(mean + min + max ~ x | coef, data=prob, type="l",
                       lty=c(1,2,2), col=c("black", "red", "red"), xlab=xlab,
                       scales=scales, ylab=ylab, as.table=as.table,
                       between=between, ...)
    }
    else {
        plot <- xyplot(pred ~ x | coef, data=prob, type="l", col="black",
                       xlab=xlab, scales=scales, ylab=ylab, as.table=as.table,
                       between=between, ...)
    }
    
    probobj <- list(est=prob, coef.labels=obj$coef.labels,
                    default.plot=plot)
    class(probobj) <- c("boolprob", class(probobj))
    probobj  
}

### ----------------------------------------------------------------------------
### Utilities
### ----------------------------------------------------------------------------

## Calculate predicted probabilities.
calc_prob <- function(newdata, conf.int, n, obj)
{    
    do_work <- function(prob_fn) {
        prob <- lapply(newdata, function(X) apply(X, 1, prob_fn))
        prob <- lapply(names(prob), function (X) {
            data.frame(pred=prob[[X]], x=newdata[[X]][,X], coef=X)
        })
        do.call(rbind, prob)    
    }
    
    if (conf.int) {
        sim_coef <- rmvnorm(n, mean=coef(obj)[,1], sigma=vcov(obj)[[1]])
        
        prob <- lapply(1:n,
                       function(X)
                       {
                           obj$model.fit[[1]]$par <- sim_coef[X,]
                           do_work(make_prob(obj))
                       })
        
        Y <- do.call(cbind, lapply(prob, function(x) x[,1]))
        q <- t(apply(Y, 1, quantile, probs=c(0.025, 0.975)))
        Y <- cbind(rowMeans(Y), q, prob[[1]][,c("x", "coef")])
        names(Y) <- c("mean", "min", "max", "x", "coef")
        prob <- Y
    }
    else {
        ## Define predicted probability function.
        prob <- do_work(make_prob(obj))
    }  
    prob    
}

## Returns newdata matrix given the suppied obj, vars, and newdata.
set_newdata <- function(obj, vars, newdata, k)
{
    ## either vars or newdata is required.
    if (all(is.null(vars), is.null(newdata)))
        stop("either 'vars' or 'newdata' must be supplied", .call=FALSE)
    
    ## Return newdata as is if specified by user.
    if (!is.null(newdata)) {
        if (!is.list(newdata) | !all(names(newdata) %in% obj$coef.labels))
            stop("newdata must be a list of data.frames and the name of list items must correspond \n   to the variable being varied", .call=FALSE)
        return(newdata)
    }
    
    ## get index values of vars if names instead of indices are supplied.
    if (is.character(vars))
        vars <- match(vars, obj$coef.labels)
    
    ## function to build sets of covariate profiles for specified variable.
    ## var.idx is variable index identifying the variable to vary; the other
    ## variables are defined in the calling environment.
    build_profile <- function(var.idx) {
        profiles <- matrix(mu, ncol=length(mu), nrow=k, byrow=TRUE)
        profiles[,var.idx] <- seq(min(mm[,var.idx]), max(mm[,var.idx]), length=k)
        profiles
    }
    
    mm  <- model.matrix(obj)
    mu  <- colMeans(mm)
    newdata <- Map(build_profile, vars)
    names(newdata) <- obj$coef.labels[vars]
    for (n in names(newdata))
        colnames(newdata[[n]]) <- obj$coef.labels
    newdata
}

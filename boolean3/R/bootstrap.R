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

##' Performs bootstrap estimates of a boolean model.
##'
##' \code{boolboot} performs bootstrap estimated of a boolean model specified
##' by \code{\link{boolprep}} and first estimated with \code{\link{boolean}}.
##' @title Bootstrap for Boolean Models
##' @param obj boolean model object as produced by \code{\link{boolprep}} and
##' first estimated with \code{\link{boolean}}.
##' @param n integer specifying the number of bootstrap estimates. Defaults to
##' \code{100}.
##' @param method string specifying the method of estimation. The specified
##' method should be one of those available from the \code{\link{optimx}} or
##' \code{\link{optim}} functions. Defaults to \code{"nlminb"}.
##' @param cluster string vector specifying hosts to use for parallel processing
##' through \code{parallel} (see \code{\link{makeCluster}}). Defaults to \code{NULL}
##' indicating no clustering.
##' @param ... additional parameters to pass on to subsequent functions.
##' @return \code{boolboot} returns a \code{boolboot} model object. This
##' object is identical to a boolean model object but with an additional
##' \code{model.boot} slot containing the results of the bootstrap. A separate
##' object type is used to help prevent the accidential loss of bootstrap
##' estimates.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @references Braumoeller, Bear F. (2003) ``Causal Complexity and the Study
##' of Politics.'' \emph{Political Analysis} 11(3): 209--233.
##' @seealso See \code{\link{boolprep}} for model setup, \code{\link{boolean}}
##' for estimation, the \code{parallel} package for details on clustering,
##' and \code{\link{optimx}} and \code{\link{optim}} for estimation methods.
##' @export
boolboot <- function(obj, n=100, method="nlminb", cluster=NULL, ...)
{    
    if (!inherits(obj, "boolean"))
        stop("boolboot requires a boolean object as outputed by boolprep", call.=FALSE)
    if (inherits(obj, "boolboot"))
        stop("model object already contains bootstrapping results. Halting to prevent data loss", call.=FALSE)        

    ## A single method has to be specified. Defaults to nlminb
    if (length(method) > 1)
        obj$method <- method[1]
    else
        obj$method <- method
    
    if (!is.null(cluster)) {
        ## Bootstrap the model with parallel cluster.
        bs_closure <- bs_mod(obj, ...)
        Results <- do_parallel(cluster, bs_closure, n, ...)
        Results <- t(do.call(rbind, Results))
    }
    else {
        ## Bootstrapping without cluster.
        bs_closure <- bs_mod(obj, ...)
        Results <- sapply(1:n, bs_closure, ...)
    }

    obj[["model.boot"]] <- list(coef=t(Results), samples=n)
    class(obj) <- c("boolboot", class(obj))
    obj
}

do_parallel <- function(cluster, fn, n, ...)
{
    Cluster <- makeCluster(cluster, type="SOCK")
    clusterSetRNGStream(Cluster)    
    Results <- clusterApplyLB(Cluster, 1:n, fn, ...)
    stopCluster(Cluster)
    Results
}

bs_sample <- function(obj, ...)
{
    ## Samples with replacement from the model.matrix of a boolean model object.
    S <- model.matrix(obj)
    Y <- obj$response
    idx <- sample(1:nrow(S), replace=TRUE)
    list(y=Y[idx], samp=S[idx,])
}

bs_obj <- function(obj, ...)
{
    ## Returns an augmented model object that includes the new bootstrapped
    ## sample.
    tmp <- bs_sample(obj)
    new.y <- tmp$y
    new.samp <- tmp$samp
    
    for (mod in names(obj$coef.idx)) {
        nam <- colnames(obj$frame[[mod]])
        obj$frame[[mod]] <- new.samp[, obj$coef.idx[[mod]]]
        colnames(obj$frame[[mod]]) <- nam
    }
    obj$response <- new.y
    obj$calc_lik <- make_lik(obj)
    obj
}

bs_mod <- function(obj, ...)
{
    ## Function used internally for bootstrapping a model. Returns a closure.
    n <- 1
    function(...) {
        cat("Bootstrap interation:", n, "\n")
        n <<- n+1
        bs_fit <- boolean(bs_obj(obj), ...)
        bs_fit$model.fit[[1]]$par
    }
}

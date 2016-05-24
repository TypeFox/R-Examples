## ----------------------------------------------------------------------------
## This file is part of boolean3
##
## Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
##
## boolean3 represents a substantial re-write of the original boolean package
## developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
## was developed under the direction of Bear Braumoeller and with support from
## The Ohio State University's College of Social and Behavioral Sciences.
##
## boolean3 and is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see <http://www.gnu.org/licenses/>.
##
## ----------------------------------------------------------------------------

##' Extract model matrix from boolean model.
##'
##' This function extracts the model matrix from the specified boolean
##' model. Note that this model wouldn't always be appropriate for estimating a
##' model since multiple intercepts are included, making the columns of the
##' matrix perfectly collinear.
##' @method model.matrix boolean
##' @title Model Matrix from boolean object
##' @param object A boolean object.
##' @param ... Arguments to pass on.
##' @return An n-by-k model matrix.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export model.matrix.boolean
model.matrix.boolean <- function(object, ...) {
  if (!("boolean" %in% class(object))) {
    stop("model.matrix requires a boolean object as outputed by boolprep", .call=FALSE)
  }

  frame <- c()
  for (f in names(object$frame)) {
    frame <- cbind(frame, object$frame[[f]])
  }

  colnames(frame) <- object$coef.labels
  frame
}

##' Calculate predicted probabilities
##'
##' Calculate predicted probabilities for each observation.
##' @method predict boolean
##' @title Predicted Probabilities
##' @param object boolean object.
##' @param newdata optionally, a data.frame containing the covariate profiles to
##' predict. If omitted, the original data used in fitting the model will be
##' used. Note that newdata must have the same structure as that returned by
##' model.matrix.
##' @param ... Additional parameters to pass on.
##' @return A vector of predicted probabilities for the covariate profiles
##' specified in newdata.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export predict.boolean
predict.boolean <- function(object, newdata=NULL, ...)
{
    if (is.null(object$model.fit) | !("boolfit" %in% class(object$model.fit)))
        stop("predict requires a boolfit object", .call=FALSE)

    if (!is.null(newdata))
        if (ncol(model.matrix(object)) != ncol(newdata))
            stop("newdata must have the same structure as the original model matrix",
                 .call=FALSE)

    obj <- object                       # Solves some namespace bug
    
    fn <- make_prob(obj)

    if (is.null(newdata))
        newdata <- model.matrix(obj)

    apply(newdata, 1, fn)    
}

##' Calculate the variance-covariance matrix for a boolean object.
##'
##' Calculate the variance-covariance matrix for a boolean object.
##' @title Variance-Covariance Matrix for Boolean Model
##' @param object boolean object.
##' @param ... Additional parameters to pass on.
##' @return A matrix, the negative of the inverse Hessian.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
vcov.boolean <- function(object, ...) {
  lapply(object$model.fit, function(x) solve(x$nhatend))
}

##' Return the log-likelihood for a boolean model.
##'
##' This function returns the log-likelihood of the estimated boolean model.
##'
##' @method logLik boolean
##' @param object boolean model object.
##' @param ... Additional parameters to be passed.
##' @return This function returns an object of class \code{\link{logLik}}.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export logLik.boolean
logLik.boolean <- function(object, ...)
{
    if(is.null(object$model.fit) | !("boolfit" %in% class(object$model.fit))) {
        stop("logLik requires a boolfit object", .call=FALSE)
    }

    ## Convert NULLs to NA; otherwise, AIC/BIC breaks.
    llik <- lapply(object$model.fit,
                   function(x)
                   {
                       if (is.null(x$value)) {
                           NA
                       }
                       else {
                           x$value
                       }
                   })

  attr(llik, "df") <- object$k
  class(llik) <- "logLik"
  llik
}

##' Extract the Number of Observations from a Fit.
##'
##' Extract the Number of Observations from a Fit.
##' @method nobs boolean
##' @title Extract the Number of Observations from a Fit.
##' @param object boolean model object
##' @param ... further arguments to be passed to other methods.
##' @return An integer.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export nobs.boolean
nobs.boolean <- function(object, ...) {
  object$N
}

##' Return the model estimates.
##'
##' This function returns a vector of coefficient estimates.
##'
##' @method coef boolean
##' @param object boolean model object.
##' @param ... Additional parameters to be passed.
##' @return A vector of the coefficient estimates.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export coef.boolean
coef.boolean <- function(object, ...) {
  if(is.null(object$model.fit) | !("boolfit" %in% class(object$model.fit))) {
      stop("coef requires a boolfit object", .call=FALSE)
  }

  ## Extract the coefficients, label them.
  sapply(object$model.fit, function(x) x$par)
}

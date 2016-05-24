# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    

# wrapper for Tps object 
"predict.Tps"<- function(object, ...){
  UseMethod("Krig")
 }

"predict.Krig" <- function(object, x = NULL, Z = NULL, 
    drop.Z = FALSE, just.fixed = FALSE, lambda = NA, df = NA, 
    model = NA, eval.correlation.model = TRUE, y = NULL, yM = NULL, 
    verbose = FALSE, ...) {
    #NOTE: most of this function is figuring out what to do!
    #
    # check that derivative is not called
    if (!is.null(list(...)$derivative)) {
        stop("For derivatives use predictDerivative")
    }
    # y is full data yM are the data collapsed to replicate means
    # if new data is not passed then copy from the object
    if (is.null(y) & is.null(yM)) {
        temp.c <- object$c
        temp.d <- object$d
    }
    # check for passed x but no Z -- this is an error
    # if there  are Z covariates in the model and drop.Z is FALSE
    ZinModel<- !is.null(object$Z)
    newX<- !is.null(x)
    missingZ<- is.null(Z)
    if( ZinModel&newX){
    if( missingZ  & !drop.Z) {
        stop("Need to specify drop.Z as TRUE or pass Z values")
    }
    }
    # default is to predict at data x's
    if (is.null(x)) {
        x <- object$x
    }
    else {
        x <- as.matrix(x)
    }
    # default is to predict at data Z's
    if (is.null(Z)) {
        Z <- object$Z
    }
    else {
        Z <- as.matrix(Z)
    }
    if (verbose) {
        print(x)
        print(Z)
    }
    # transformations of x values used in Krig
    xc <- object$transform$x.center
    xs <- object$transform$x.scale
    x <- scale(x, xc, xs)
    # NOTE knots are already scaled in Krig object and are used
    # in transformed scale.
    #  i.e.   knots <- scale( object$knots, xc, xs)
    #
    # figure out if the coefficients for the surface needto be recomputed.
    find.coef <- (!is.null(y) | !is.null(yM) | !is.na(lambda) | 
        !is.na(df) | !is.na(model[1]))
    if (verbose) {
        cat("find.coef", find.coef, fill = TRUE)
    }
    #   convert effective degrees of freedom to equivalent lambda
    if (!is.na(df)) {
        lambda <- Krig.df.to.lambda(df, object$matrices$D)
    }
    if (!is.na(model)) {
        lambda <- model[1]
    }
    if (is.na(lambda)) 
        lambda <- object$lambda
    #
    # if the coefficients need to be recomputed  do it.
    if (find.coef) {
        if (verbose) {
            cat("new coefs found", fill = TRUE)
        }
        object3 <- Krig.coef(object, lambda = lambda, y = y, 
            yM = yM)
        temp.d <- object3$d
        temp.c <- object3$c
    }
    if (verbose) {
        cat(" d coefs", fill = TRUE)
        print(temp.d)
        cat("c coefs", fill = TRUE)
        print(temp.c)
    }
    
    # this is the fixed part of predictor
    #
    Tmatrix <- do.call(object$null.function.name, c(object$null.args, 
        list(x = x, Z = Z, drop.Z = drop.Z)))
    if (drop.Z) {
        temp <- Tmatrix %*% temp.d[object$ind.drift]
    }
    else {
        temp <- Tmatrix %*% temp.d
    }
    # add in spatial piece
    if (!just.fixed) {
        #
        # Now find sum of covariance functions times coefficients
        # Note that the multiplication of the cross covariance matrix
        # by the coefficients is done implicitly in the covariance function
        #
        # The covariance function is
        # evaluated by using its name, the do.call function, and any
        # additional arguments.
        #
        temp <- temp + do.call(object$cov.function.name, c(object$args, 
            list(x1 = x, x2 = object$knots, C = temp.c)))
    }
    #
    # transform back into raw scale if this is a correlation model.
    # if y's are in the scale of correlations
    # if so scale by sd and add back in mean
    correlation.model <- (object$correlation.model & eval.correlation.model)
    if (correlation.model) {
        if (!is.na(object$sd.obj[1])) {
            temp <- temp * predict(object$sd.obj, x)
        }
        if (!is.na(object$mean.obj[1])) {
            temp <- temp + predict(object$mean.obj, x)
        }
    }
    return(temp)
}

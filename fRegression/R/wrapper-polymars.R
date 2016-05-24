
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# INTERFACE:            FROM POLSPLINE - POLYMARS DESCRIPTION:
#  .polymarsFormula      Polymars regress from package polspline
#  .polymars.default     Default wrapper for polymars()
#  .predict.polymars     Formula wrapper for polymars()
#  .predict.polymars     Predict from a polymars model
################################################################################


# Note:
    # Introduce no .polymars = function() UseMethod()
    #   this fails regFit(..., use = "polymars)


# ------------------------------------------------------------------------------


.polymarsFormula <- 
    function(formula, data, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Extract Model Data:
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)

    # Rempove Intercept from x if exists ...
    M <- which(colnames(x) == "(Intercept)")
    if (length(M) > 0) X <- x[ ,-M]

    # Fit:
    fit <- .polymarsDefault(responses = y, predictors = X, ...)

    # Add to fit:
    #   ... '$coef' keeps model
    fit$model <- mf
    fit$terms <- mt

    # Class:
    class(fit) <- "polymars"

    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


.polymarsDefault <-
    function(responses, predictors, maxsize, gcv = 4, additive = FALSE,
    startmodel, weights, no.interact, knots, knot.space = 3, ts.resp, ts.pred,
    ts.weights, classify, factors, tolerance = 1e-06, verbose = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #  responses - a vector (or matrix) of responses. (Can be a a vector of
    #       characters for classification)
    #  predictors - a matrix of predictors with same number of cases as
    #       response. Columns are predictors.

    # Optional Arguments:
    #  maxsize - maximum number of basis function the model can contain
    #  gcv - parameter for overall best model seletion
    #  additive - boolean, is the model to be additive
    #  startmodel - either a matrix (m*4 or m*5) or a polymars object from
    #       a previous call to polymars
    #       an initial model the procedure should start with in model
    #       selection
    #  weights - a vector of length equal to the number of cases
    #  no.interact - a 2*l matrix of columns numbers of the predictor
    #       matrix (each row pair cannot have interaction terms)
    #  knots - a vector specifying many knots per predictor are
    #       wanted (with -1 for categorical variables)
    #       ncol(predictors)==length(knots), or a matrix with
    #       ncol(predictors) == ncol(knots) with actual knot
    #       specified and filled out with NA's.
    #       Can also be a single number - "knots" number of knots
    #       per predictor
    #  knot.space - minimum number of order statistics between knots
    #  ts.resp - testset reponses, same format as responses
    #  ts.pred - testset predictors, same format as predictors
    #  ts.weights - testset weights, same format as weights
    #  classify - whether classification is to be done, set = TRUE if the
    #       response vector is integer, if
    #       if character classify is automatically true
    #  factors - a vector of column numbers of the predictor matrix of
    #       categorical variables
    #  tolerance - a numerical parameter which may need to be made smaller
    #       if the program crashes store the call to the polymars
    #       function

    # FUNCTION:

    # require(polspline)

    print(head(responses))
    print(head(predictors))

    # Fit:
    .Call <- match.call()
    .Call[[1]] <- quote(polspline::polymars)
    ans <- eval(.Call, parent.frame())

    # Add Coefficients Parameters:
    ans$coef <- ans$model
    ans$parameters <- ans$coef
    ans$fitted.values <- ans$fitted

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.predict.polymars <-
    function(object, newdata, se.fit = FALSE, type = "response", ...)
{
    # Note:
    #   newdata is a predictor data.frame, if missing the fitted
    #   vector will be returned.

    # Example:
    #   x=LM3(); object1 = regFit(Y ~ X1+X2+X3, data = x, use = "polymars")@fit
    #   .predict.polymars(object, newdata = x[, -1])

    # FUNCTION:

    # Restore Object Model:
    object$model <- object$coef
    class(object) <- "polymars"

    # Polymars requires 1-column matrices:
    object$residuals <- matrix(object$residuals)
    object$fitted <- matrix(object$fitted)

    # Here, object is expected to be the slot @fit of an object of class 'fREG'
    if (missing(newdata)) {
        y <- as.vector(object$fitted)
    } else {
        tt <- object$terms
        Terms <- delete.response(tt)
        modelFrame <- model.frame(Terms, newdata)
        X <- model.matrix(Terms, modelFrame)[, -1]
        Y <- polspline::predict.polymars(object, x = X, ...)
    }

    # Add optionally standard errors - NA's not available yet ...
    if (se.fit) Y <- list(fit = Y, se.fit = NA*Y)

    # Return Value:
    Y
}


################################################################################



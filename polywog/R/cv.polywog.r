## Declare `d` as a global variable so R CMD check doesn't complain about the
## foreach loop
if (getRversion() >= "2.15.1")
    utils::globalVariables("d")

##' k-fold cross-validation to select the polynomial degree and penalization
##' factor for a \code{\link{polywog}} model.
##'
##' When fitting with \code{method = "scad"}, different fold assignments are
##' used for each polynomial degree specified, because \code{\link{cv.ncvreg}}
##' does not allow for custom fold assignments.  This may affect the accuracy
##' of the estimated cross-validation error for each degree.  When
##' \code{method = "scad"}, the calls to \code{\link{polywog}} made by
##' \code{cv.polywog} will issue warnings that the \code{foldid} argument is
##' being ignored.
##' @title Cross-validation of polynomial degree and penalization factor
##' @param formula model formula specifying the response and input variables.
##' @param ... other arguments to be passed to \code{\link{polywog}}.  Arguments
##' related to the bootstrap will be ignored, as bootstrapping must be
##' performed separately.
##' @param degrees.cv vector of polynomial degrees to examine via
##' cross-validation.
##' @param nfolds number of folds to use in cross-validation.
##' @param model logical: whether to include the model frame in the
##' \code{"polywog"} object included in the output.
##' @param X logical: whether to include the raw model matrix (i.e., the
##' matrix of input variables prior to taking their polynomial expansion) in
##' the \code{"polywog"} object included in the output.
##' @param y logical: whether to include the response variable in the
##' \code{"polywog"} object included in the output.
##' @return An object of class \code{"cv.polywog"}, a list containing:
##' \describe{
##'   \item{\code{results}}{A table of each degree tested, the optimal
##' penalization factor \eqn{\lambda} for that degree, and its
##' cross-validation error.}
##'   \item{\code{degree.min}}{The polynomial degree giving the lowest
##' cross-validation error.}
##'   \item{\code{polywog.fit}}{A \code{\link{polywog}} model, fit at the
##' polynomial degree giving the lowest cross-validation error.}
##' }
##'
##' Because the returned object contains the fitted polywog model for the
##' optimal degree, no additional runs of \code{\link{polywog}} are necessary
##' to estimate coefficients or the penalization factor \eqn{\lambda}.
##' However, bootstrap results must be obtained by running
##' \code{\link{bootPolywog}} on the \code{"polywog.fit"} element of the
##' returned object, as in the examples below.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
##' @example inst/examples/cv.polywog.r
##' @import foreach
cv.polywog <- function(formula,
                       ...,
                       degrees.cv = 1:3,
                       nfolds = 10,
                       model = TRUE,
                       X = FALSE,
                       y = FALSE)
{
    cl <- match.call()

    ## Assign cross-validation folds
    ##
    ## Unfortunately, this requires knowing the number of observations in the
    ## dataset in advance.  We could do a trial run of polywog() with high
    ## tolerance/low maxit to get that, but that still entails taking the QR
    ## decomposition of the model matrix to check for singularities, which can
    ## be costly with large N.  So instead we form (and immediately discard)
    ## the model frame.
    if (nfolds < 2)
        stop("'nfolds' must be at least 2")
    formula <- as.Formula(formula)
    nobs <- match(c("data", "subset", "na.action"), names(cl), 0L)
    nobs <- cl[c(1L, nobs)]
    nobs$formula <- formula
    nobs[[1]] <- quote(stats::model.frame)
    nobs <- nrow(eval(nobs, parent.frame()))
    foldid <- sample(seq_len(nfolds), size = nobs, replace = TRUE)

    ## Run polywog() on each of the specified degrees
    polyFits <- foreach(d = degrees.cv) %do% {
        ## Construct a call to polywog() using the specified degree
        ##
        ## It would be possible to instead do this via
        ##
        ##   dots <- list(...)
        ##   dots$formula <- formula
        ##   dots$degree <- d
        ##   ## etc
        ##   do.call("polywog", dots)
        ##
        ## but then the "call" element of the fitted model would contain the
        ## evaluated arguments instead of their names, so printing it or
        ## calling summary() on it would end up printing the entire dataset
        ## used to fit the model
        polyCall <- cl
        polyCall$degree <- d
        polyCall$foldid <- foldid
        polyCall$degrees.cv <- polyCall$nfolds <- NULL
        polyCall$boot <- 0
        polyCall$model <- polyCall$X <- polyCall$y <- FALSE
        polyCall[[1]] <- quote(polywog)

        eval(polyCall, parent.frame())
    }

    ## Construct table of optimal lambda and minimal MSE by degree
    results <- cbind(degree = degrees.cv,
                     lambda.min = sapply(polyFits, "[[", "lambda"),
                     cv.err = sapply(polyFits, function(x) x$lambda.cv$errorMin))

    ## Find the optimal degree
    indMin <- which.min(results[, "cv.err"])
    degreeMin <- degrees.cv[indMin]

    ## Construct the fitted model object corresponding to the optimal degree
    fitMin <- polyFits[[indMin]]
    fitMin$call$foldid <- NULL
    if (model) {
        fitMin$model <- model.frame(fitMin)
        fitMin$call$model <- TRUE
    }
    if (X) {
        fitMin$X <- model.matrix(fitMin)
        fitMin$call$X <- TRUE
    }
    if (y) {
        fitMin$y <- model.response(model.frame(fitMin))
        fitMin$call$y <- TRUE
    }

    structure(list(results = results,
                   degree.min = degreeMin,
                   polywog.fit = fitMin),
              class = "cv.polywog")
}

##' @useDynLib polywog
##' @importFrom Rcpp sourceCpp
NULL

##' Bootstrapped basis regression with oracle model selection
##'
##' A package for flexible functional form estimation via bootstrapped basis
##' regression with oracle model selection.  This version of the software should
##' be considered \strong{in beta}.  For bug reports and feature requests,
##' please email Brenton Kenkel (\email{brenton.kenkel@@gmail.com}) or file an
##' issue at \url{https://github.com/brentonk/polywog-package/issues}.
##' @name polywog-package
##' @docType package
##' @section Acknowledgements: We are grateful to Tyson Chatangier for many
##' helpful suggestions about an earlier version of the package.  We also thank
##' the Wallis Institute of Political Economy and the Theory and Statistics
##' Research Lab at the University of Rochester for financial support during the
##' writing of the package.
##' @references
##' Brenton Kenkel and Curtis S. Signorino.  2012.  "A Method for Flexible
##' Functional Form Estimation: Bootstrapped Basis Regression with Variable
##' Selection."  Typescript, University of Rochester.
NULL

##' Polynomial regression with oracle variable selection
##'
##' Fits a regression model using a polynomial basis expansion of the input
##' variables, with penalization via the adaptive LASSO or SCAD to provide
##' oracle variable selection.
##'
##' The design matrix for the regression is a polynomial basis expansion of the
##' matrix of raw input variables.  This includes all powers and interactions of
##' the input variables up to the specified \code{degree}.  For example, the
##' following terms will be included in \code{polywog(y ~ x1 + x2, degree = 3,
##' ...)}:
##' \itemize{
##'   \item terms of degree 0: intercept
##'   \item terms of degree 1: \code{x1}, \code{x2}
##'   \item terms of degree 2: \code{x1^2}, \code{x2^2}, \code{x1*x2}
##'   \item terms of degree 3: \code{x1^3}, \code{x2^3}, \code{x1*x2^2},
##' \code{x1^2*x2}
##' }
##' To exclude certain terms from the basis expansion, use a model formula like
##' \code{y ~ x1 + x2 | z1 + z2}.  Only the degree 1 terms of \code{z1} and
##' \code{z2} will be included.
##'
##' It is possible that the "raw" basis expansion will be rank-deficient, such
##' as if there are binary input variables (in which case \eqn{x_i = x_i^n} for
##' all \eqn{n > 0}).  The procedure detects collinearity via \code{\link{qr}} and
##' removes extraneous columns before fitting.
##'
##' For both the adaptive LASSO and SCAD, the penalization factor \eqn{\lambda}
##' is chosen by k-fold cross-validation.  The selected value minimizes the
##' average mean squared error of out-of-sample fits.  (To select both
##' \eqn{\lambda} and the polynomial degree simultaneously via cross-validation,
##' see \code{\link{cv.polywog}}.)
##'
##' The cross-validation process may be run in parallel via
##' \code{\link{foreach}} by registering an appropriate backend and specifying
##' \code{.parallel = TRUE}.  The appropriate backend is system-specific; see
##' \code{\link{foreach}} for information on selecting and registering a
##' backend.  The bootstrap iterations may also be run in parallel by
##' specifying \code{control.boot = control.bp(.parallel = TRUE)}.
##' @param formula model formula specifying the response and input
##' variables.  See "Details" for more information.
##' @param data a data frame, list or environment containing the variables
##' specified in the model formula.
##' @param subset an optional vector specifying a subset of observations to be
##' used in fitting.
##' @param weights an optional vector specifying weights for each observation to
##' be used in fitting.
##' @param na.action a function specifying what to do with observations
##' containing \code{NA}s (default \code{\link{na.omit}}).
##' @param degree integer specifying the degree of the polynomial expansion of
##' the input variables.
##' @param family \code{"gaussian"} (default) or \code{"binomial"} for logistic
##' regression (binary response only).
##' @param method variable selection method: \code{"alasso"} (default) for
##' adaptive LASSO or \code{"scad"} for SCAD.  You can also select \code{method
##' = "none"} to return the model matrix and other information without fitting.
##' @param penwt.method estimator for obtaining first-stage estimates in
##' logistic models when \code{method = "alasso"}: \code{"lm"} (default) for a
##' linear probability model, \code{"glm"} for logistic regression.
##' @param unpenalized names of model terms to be exempt from the adaptive
##' penalty (only available when \code{method = "alasso"}).
##' @param .parallel logical: whether to perform k-fold cross-validation in
##' parallel (only available when \code{method = "alasso"}).  See "Details"
##' below for more information on parallel computation.
##' @param boot number of bootstrap iterations (0 for no bootstrapping).
##' @param control.boot list of arguments to be passed to
##' \code{\link{bootPolywog}} when bootstrapping; see \code{\link{control.bp}}.
##' @param lambda a vector of values from which the penalty factor is to be
##' selected via k-fold cross-validation.  \code{lambda} is left unspecified
##' by default, in which case a sequence of values is generated automatically,
##' controlled by the \code{nlambda} and \code{lambda.min.ratio} arguments.
##' Naturally, k-fold cross-validation is skipped if \code{lambda} contains
##' exactly one value.
##' @param nlambda number of values of the penalty factor to examine via
##' cross-validation if \code{lambda} is not specified in advance; see
##' "Details".
##' @param lambda.min.ratio ratio of the lowest value to the highest in the
##' generated sequence of values of the penalty factor if \code{lambda} is
##' not specified; see "Details".
##' @param nfolds number of folds to use in cross-validation to select the
##' penalization factor.
##' @param foldid optional vector manually assigning fold numbers to each
##' observation used for fitting (only available when \code{method =
##' "alasso"}).
##' @param thresh convergence threshold, passed as the \code{thresh} argument
##' to \code{\link{glmnet}} when \code{method = "alasso"} and as the
##' \code{eps} argument to \code{\link{ncvreg}} when \code{method = "scad"}.
##' @param maxit maximum number of iterations to allow in adaptive LASSO or
##' SCAD fitting.
##' @param model logical: whether to include the model frame in the returned
##' object.
##' @param X logical: whether to include the raw design matrix (i.e., the
##' matrix of input variables prior to taking their polynomial expansion) in
##' the returned object.
##' @param y logical: whether to include the response variable in the returned
##' object.
##' @return An object of class \code{"polywog"}, a list containing: \describe{
##'   \item{\code{coefficients}}{the estimated coefficients.}
##'   \item{\code{lambda}}{value of the penalty factor \eqn{\lambda} used to
##' fit the final model.}
##'   \item{\code{lambda.cv}}{a list containing the results of the
##' cross-validation procedure used to select the penalty factor: \describe{
##'     \item{\code{lambda}}{values of the penalty factor tested in
##' cross-validation.}
##'     \item{\code{cvError}}{out-of-fold prediction error corresponding to
##' each value of \code{lambda}.}
##'     \item{\code{lambdaMin}}{value of \code{lambda} with the minimal
##' cross-validation error.}
##'     \item{\code{errorMin}}{minimized value of the cross-validation error.}
##'   }}
##'   \item{\code{fitted.values}}{the fitted mean values for each observation
##' used in fitting.}
##'   \item{\code{lmcoef}}{coefficients from an unpenalized least-squares
##' regression of the response variable on the polynomial expansion of the
##' input variables.}
##'   \item{\code{penwt}}{adaptive weight given to each term in the LASSO
##' penalty (\code{NULL} for models fit via SCAD).}
##'   \item{\code{formula}}{model formula, as a \code{\link{Formula}} object.}
##'   \item{\code{degree}}{degree of the polynomial basis expansion.}
##'   \item{\code{family}}{model family, \code{"gaussian"} or
##' \code{"binomial"}.}
##'   \item{\code{weights}}{observation weights if specified.}
##'   \item{\code{method}}{the specified regularization method.}
##'   \item{\code{penwt.method}}{the specified method for calculating
##' the adaptive LASSO weights (\code{NULL} for models fit via SCAD).}
##'   \item{\code{unpenalized}}{logical vector indicating which terms were not
##' included in the LASSO penalty.}
##'   \item{\code{thresh}}{convergence threshold used in fitting.}
##'   \item{\code{maxit}}{iteration limit used in fitting.}
##'   \item{\code{terms}}{the \code{\link{terms}} object used to construct the
##' model frame.}
##'   \item{\code{polyTerms}}{a matrix indicating the power of each raw input
##' term (columns) in each term of the polynomial expansion used in fitting
##' (rows).}
##'   \item{\code{nobs}}{the number of observations used to fit the model.}
##'   \item{\code{na.action}}{information on how \code{NA} values in the input
##' data were handled.}
##'   \item{\code{xlevels}}{levels of factor variables used in fitting.}
##'   \item{\code{varNames}}{names of the raw input variables included in the
##' model formula.}
##'   \item{\code{call}}{the original function call.}
##'   \item{\code{model}}{if \code{model = TRUE}, the model frame used in
##' fitting; otherwise \code{NULL}.}
##'   \item{\code{X}}{if \code{X = TRUE}, the raw model matrix (i.e., prior to
##' taking the polynomial expansion); otherwise \code{NULL}.  For calculating
##' the expanded model matrix, see \code{\link{model.matrix.polywog}}.}
##'   \item{\code{y}}{if \code{y = TRUE}, the response variable used in
##' fitting; otherwise \code{NULL}.}
##'   \item{\code{boot.matrix}}{if \code{boot > 0}, a sparse matrix of class
##' \code{"\linkS4class{dgCMatrix}"} where each column is the estimate from a
##' bootstrap replicate.  See \code{\link{bootPolywog}} for more information
##' on bootstrapping.}
##' }
##' @seealso To estimate variation via the bootstrap, see
##' \code{\link{bootPolywog}}.  To generate fitted values, see
##' \code{\link{predVals}} (and the underlying method
##' \code{\link{predict.polywog}}).  For plots, see \code{\link{plot.polywog}}.
##' The polynomial degree may be selected via cross-validation using
##' \code{\link{cv.polywog}}.
##'
##' Adaptive LASSO estimates are provided via \code{\link{glmnet}} and
##' \code{\link{cv.glmnet}} from the \pkg{glmnet} package.  SCAD estimates are
##' via \code{\link{ncvreg}} and \code{\link{cv.ncvreg}} in the \pkg{ncvreg}
##' package.
##' @references Brenton Kenkel and Curtis S. Signorino.  2012.  "A Method for
##' Flexible Functional Form Estimation: Bootstrapped Basis Regression with
##' Variable Selection."  Typescript, University of Rochester.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
##' @example inst/examples/polywog.r
##' @import Formula
polywog <- function(formula,
                    data,
                    subset,
                    weights,
                    na.action,
                    degree = 3,
                    family = c("gaussian", "binomial"),
                    method = c("alasso", "scad"),
                    penwt.method = c("lm", "glm"),
                    unpenalized = character(0),
                    .parallel = FALSE,
                    boot = 0,
                    control.boot = control.bp(.parallel = .parallel),
                    lambda = NULL,
                    nlambda = 100,
                    lambda.min.ratio = 1e-4,
                    nfolds = 10,
                    foldid = NULL,
                    thresh = ifelse(method == "alasso", 1e-7, 0.001),
                    maxit = ifelse(method == "alasso", 1e5, 5000),
                    model = TRUE,
                    X = FALSE,
                    y = FALSE)
{
    cl <- match.call()
    family <- match.arg(family)
    method <- match.arg(method)
    penwt.method <- match.arg(penwt.method)
    ret.X <- X
    ret.y <- y

    ## Check for bad argument combinations
    if (method == "scad" && penwt.method != "lm")
        warning("Argument 'penwt.method' is ignored when method = \"scad\"")
    if (method == "scad" && !is.null(foldid))
        warning("Argument 'foldid' is ignored when method = \"scad\"")
    if (method == "scad" && .parallel)
        warning("Cross-validation cannot be performed in parallel when method =\"scad\"")
    if (method == "scad" && length(unpenalized) > 0)
        warning("Argument 'unpenalized' is ignored when method = \"scad\"")

    ## Assemble the model frame the usual way
    formula <- as.Formula(formula)
    mf <- match(c("data", "subset", "weights", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formula
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    terms <- attr(mf, "terms")

    ## Compute the (raw) model matrix, response variable, and polynomial term
    ## matrix
    X <- makeX(formula, mf)
    y <- model.response(mf)
    if (family == "binomial" && !all(y %in% 0:1))
        stop("Response must be binary when family = \"binomial\"")
    nobs <- nrow(X)
    polyTerms <- makePolyTerms(degree = degree,
                               k_expand = attr(X, "k_expand"),
                               k_lin = attr(X, "k_lin"),
                               binary_cols = attr(X, "binary_cols"),
                               names. = colnames(X))

    ## Extract weights if any (same code as in 'lm')
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w))
        stop("'weights' must be a numeric vector")
    nowt <- is.null(w)
    if (!nowt && method == "scad")
        stop("Weights not allowed with method = \"scad\"")
    if (nowt)
        w <- rep(1, nobs)
    if (any(w < 0))
        stop("Negative weights not allowed")

    ## Store variable names (for use in 'print.polywog')
    varNames <- extractVarNames(formula, mf)

    ## Check for unpenalized columns
    unpenalized <- makeUnpenalized(unpenalized, polyTerms)

    ## Compute linear model coefficients and remove perfectly collinear terms
    ## from the polynomial expansion and the list of unpenalized terms
    lmcoef <- computeLinearCoef(X, y, polyTerms, w)
    pivot <- attr(lmcoef, "pivot")
    polyTerms <- polyTerms[pivot, , drop = FALSE]
    unpenalized <- unpenalized[pivot]

    ## Compute penalty weights (NULL when using SCAD)
    penwt <- computePenaltyWeights(X = X,
                                   y = y,
                                   weights = w,
                                   polyTerms = polyTerms,
                                   lmcoef = lmcoef,
                                   method = method,
                                   penwt.method = penwt.method,
                                   family = family,
                                   unpenalized = unpenalized)

    ## Compute model fit
    fitPolywog <- switch(method,
                         alasso = fitAdaptiveLASSO,
                         scad = fitSCAD)
    fit <- fitPolywog(X = X,
                      y = y,
                      weights = w,
                      polyTerms = polyTerms,
                      family = family,
                      penwt = penwt,
                      lambda = lambda,
                      nlambda = nlambda,
                      lambda.min.ratio = lambda.min.ratio,
                      nfolds = nfolds,
                      foldid = foldid,
                      thresh = thresh,
                      maxit = maxit,
                      .parallel = .parallel)

    ## Compute in-sample fitted values
    fittedVals <- computePredict(X = X,
                                 poly_terms = polyTerms,
                                 coef = list(main = fit$coef),
                                 forPredVals = FALSE,
                                 interval = FALSE,
                                 bag = FALSE,
                                 level = 0,
                                 transform = family == "binomial")$fit

    ## Assemble object to return.  Ordering schema is as follows:
    ##   (1) Final estimates (including penalization factor)
    ##   (2) Auxiliary quantities computed from those estimates
    ##   (3) First-stage estimates (including penalty weights)
    ##   (4) Information used in fitting (degree, observation weights, penalty
    ##   weight method, etc.)
    ##   (5) Other auxiliary information to be used for bootstrapping,
    ##   prediction, or marginal effects
    ##   (6) Other auxiliary information used for print/summary methods
    ##   (7) Optional components (model, X, y)
    ans <- list(coefficients = fit$coef,
                lambda = fit$lambda,
                lambda.cv = fit$lambda.cv,
                ## (2)
                fitted.values = fittedVals,
                ## (3)
                lmcoef = lmcoef,
                penwt = penwt,
                ## (4)
                formula = formula,
                degree = degree,
                family = family,
                weights = if (nowt) NULL else w,
                method = method,
                penwt.method = if (method == "alasso") penwt.method,
                unpenalized = if (method == "alasso") unpenalized,
                thresh = thresh,
                maxit = maxit,
                ## (5)
                terms = terms,
                polyTerms = polyTerms,
                nobs = nobs,
                na.action = attr(mf, "na.action"),
                xlevels = .getXlevels(terms, mf),
                ## (6)
                varNames = varNames,
                call = cl,
                ## (7)
                model = if (model) mf,
                X = if (ret.X) X,
                y = if (ret.y) y)
    class(ans) <- "polywog"

    ## Bootstrapping, if requested
    if (boot > 0) {
        ans$boot.matrix <- do.call(bootPolywog,
                                   c(control.boot,
                                     list(model = ans,
                                          nboot = boot,
                                          .matrixOnly = TRUE)))
    }

    return(ans)
}

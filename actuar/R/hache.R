### ===== actuar: An R Package for Actuarial Science =====
###
### Credibility in the regression case using the Hachemeister (1975)
### model with possibly an adjustment to put the intercept at the
### barycenter of time (see Buhlmann & Gisler, 2005).
###
### AUTHORS: Xavier Milhaud, Tommy Ouellet, Vincent Goulet
### <vincent.goulet@act.ulaval.ca>

hache <- function(ratios, weights, formula, data, adj.intercept = FALSE,
                  method = c("unbiased", "iterative"),
                  tol = sqrt(.Machine$double.eps),
                  maxit = 100, echo = FALSE)
{
    Call <- match.call()

    ## If weights are not specified, use equal weights as in
    ## Buhlmann's model.
    if (missing(weights))
    {
        if (any(is.na(ratios)))
            stop("missing ratios not allowed when weights are not supplied")
        weights <- array(1, dim(ratios))
    }

    ## Check other bad arguments.
    if (NCOL(ratios) < 2)
        stop("there must be at least one node with more than one period of experience")
    if (NROW(ratios) < 2)
        stop("there must be more than one node")
    if (!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in 'weights' and in 'ratios'")
    if (all(!weights, na.rm = TRUE))
        stop("no available data to fit model")

    ## Build the design matrix
    mf <- model.frame(formula, data, drop.unused.levels = TRUE)
    mt <- attr(mf, "terms")
    xreg <- model.matrix(mt, mf)

    ## Do computations in auxiliary functions.
    res <-
        if (adj.intercept)
            hache.barycenter(ratios, weights, xreg,
                             method = match.arg(method),
                             tol = tol, maxit = maxit, echo = echo)
        else
            hache.origin(ratios, weights, xreg,
                         tol = tol, maxit = maxit, echo = echo)

    ## Add the terms object to the result for use in predict.hache()
    ## [and thus predict.lm()].
    res$terms <- mt

    ## Results
    attr(res, "class") <- "hache"
    attr(res, "model") <- "regression"
    res
}

predict.hache <- function(object, levels = NULL, newdata, ...)
{
    ## If model was fitted at the barycenter of time (there is a
    ## transition matrix in the object), then also convert the
    ## regression coefficients in the base of the (original) design
    ## matrix.
    if (!is.null(R <- object$transition))
    {
        for (i in seq_along(object$adj.models))
        {
            b <- coefficients(object$adj.models[[i]])
            object$adj.models[[i]]$coefficients <- solve(R, b)
        }
    }

    ## Prediction (credibility premiums) using predict.lm() on each of
    ## the adjusted individual models. This first requires to add a
    ## 'terms' component to each adjusted model.
    f <- function(z, ...)
    {
        z$terms <- object$terms
        class(z) <- "lm"                # to keep predict.lm() quiet
        unname(predict.lm(z, ...))
    }

    structure(sapply(object$adj.models, f, newdata = newdata), ...)
}

print.hache <- function(x, ...)
    print.default(x)

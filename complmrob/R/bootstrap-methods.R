#' Bootstrap statistics functions
#'
#' Functions to calculate the coefficient(s) of the robust linear regression model
#' from a bootstrapped sample
#'
#' Different approaches for bootstrapping have been implemented. The default "fast and robust bootstrap"
#' (FRB) proposed by M. Salibian-Barrera, et al. (2002), implemented with \code{bootStatFast} is the
#' fastest and most resistant to outliers, while the other two \code{bootStatResiduals} and \code{bootStatCases}
#' are standard bootstrap methods, where the residuals resp. the cases are resampled and the model is
#' fit to this data.
#'
#' @param origData the original data set.
#' @param residData the original data set with the columns fit, resid and the predictor variables instead of
#'      the response variable.
#' @param model The lmrob model
#' @param inds the resampled indices.
#' @param coefind the index of the coefficient to extract.
#' @param control the control object as returned by \code{bootStatFastControl}.
#' @param formula the formula to fit the model
#' @param intercept if the model includes an intercept term.
#' @param maxTries the maximum number of tries to increase the maxit control arguments for the S estimator.
#'
#' @references M. Salibian-Barrera, S. Aelst, and G. Willems. Fast and robust bootstrap. Statistical Methods and Applications, 17(1):41-71, 2008.
#' @seealso \code{\link{bootcoefs}}
#' @name bootStat-methods
NULL

#' @rdname bootStat-methods
#' @import robustbase
bootStatResiduals <- function(residData, inds, coefind, intercept = TRUE, maxTries = 4L) {
    formula <- fit + resid[inds] ~ .;

    if(intercept == FALSE) {
        formula <- fit + resid[inds] ~ . - 1;
    }

    suppressWarnings(m <- robustbase::lmrob(formula, data = residData));

    # Ensure convergence!
    itcount <- 0L;
    while(!m$converged && itcount < maxTries) {
        itcount <- itcount + 1L;
        suppressWarnings(m <- update(m, maxit.scale = m$control$maxit.scale + 200, max.it = m$control$max.it + 50))
    }
    if(itcount == maxTries) {
        return(NA_real_);
    }
    return(coef(m)[coefind])
};

#' @import robustbase
#' @rdname bootStat-methods
bootStatCases <- function(origData, inds, coefind, formula, maxTries = 4L) {
    suppressWarnings(m <- robustbase::lmrob(formula, data = origData[inds, ]));

    # Ensure convergence!
    itcount <- 0L;
    while(!m$converged && itcount < maxTries) {
        itcount <- itcount + 1L;
        suppressWarnings(m <- update(m, maxit.scale = m$control$maxit.scale + 200, max.it = m$control$max.it + 50))
    }
    if(itcount == maxTries) {
        return(NA_real_);
    }
    return(coef(m)[coefind])
}

#' @import robustbase
#' @rdname bootStat-methods
bootStatFastControl <- function(model) {
    X <- model.matrix(model);
    n <- attr(X, "dim")[1];
    p <- attr(X, "dim")[2];
    vfactorinv <- (n - p) * model$init.S$control$bb;

    scaledRes <- model$residuals / model$scale;
    scaledResS <- model$init.S$residuals / model$scale;

    w <- model$rweights / model$scale;
    v <- (model$scale / vfactorinv) *
        robustbase::Mchi(scaledResS, deriv = 0, cc = model$init.S$control$tuning.chi, psi = model$init.S$control$psi);

    wp2 <- robustbase::Mpsi(scaledRes, deriv = 1, cc = model$control$tuning.psi, psi = model$control$psi);

    sf <- crossprod(X, diag(wp2)) %*% X;
    sfinv <- NULL;

    tryCatch({
        sfinv <- solve(sf);
    }, error = function(e) {
        warning("The data is (almost) singular, results will not be reliable.");
        svddecomp <- svd(sf);
        sfinv <<- tcrossprod(svddecomp$v %*% diag(1/svddecomp$d), svddecomp$u);
    });
    M <- model$scale * sfinv %*% crossprod(X, diag(w)) %*% X

    chideriv <- robustbase::Mchi(scaledResS, deriv = 1, cc = model$init.S$control$tuning.chi, psi = model$init.S$control$psi);

    a <- (chideriv %*% scaledResS)[1, 1, drop = TRUE];

    # the correction vector d is not divided by the scale in the paper, but in the
    # reference code from Matias Salibian-Barrera
    d <- (sfinv %*% crossprod(X, wp2 * model$residuals)) * vfactorinv / (a * model$scale);

    ## The .lm.fit method is only available in R 3.1.0 upwards,
    ## so write a wrapper that can be used in the bootStatFast method
    ## that doesn't do any checks

    lmwfit <- function(x, y, weights) {};

    if(getRversion() < numeric_version("3.1.0")) {
        lmwfit <- function(x, y, weights) {
            return(lm.wfit(x, y, w = weights)$coefficients);
        };
    } else {
        # the fast version w/o any checks is available
        w <- sqrt(w);
        lmwfit <- function(x, y, weights) {
            return(.lm.fit(x * weights, y * weights)$coefficients);
        };
    }

    ret <- list(
        M = M,
        d = d,
        coefficients = coef(model),
        v = v,
        wts = w,
        scale = model$scale,
        residualsS = model$init.S$residuals,
        terms = terms(model),
        lmwfit = lmwfit
    )

    class(ret) <- "bootStatFastControl";
    return(ret);
}

#' @import robustbase
#' @rdname bootStat-methods
bootStatFast <- function(origData, inds, control, coefind) {
    bsResidS <- control$residualsS[inds];

    mf <- model.frame(control$terms, origData[inds, , drop = FALSE]);
    X <- model.matrix(control$terms, mf);
    y <- model.response(mf, "numeric");

    wts <- control$wts[inds];

    bsCoef <- control$lmwfit(X, y, wts);
    bsScale <- sum(control$v[inds]);

    bootCoefs <- control$coefficients + control$M %*% (bsCoef - control$coefficients) + control$d * (bsScale - control$scale);

    return(bootCoefs[coefind])
}

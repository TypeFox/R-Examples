adaptive.weights <- function(x, y, nu=1, weight.method=c("multivariate", "univariate")) {

    weight.method <- match.arg(weight.method)

    nobs <- nrow(x)
    if (nobs+1 <= ncol(x)) {
        warning("using univariate weight method since p>n")
        weight.method <- "univariate"
    }

    # Get OLS fits
    weights <- switch(weight.method,
                      "univariate" = apply(x, 2, function(xi){lm.fit(x=cbind(1,xi), y=y)$coefficients[2]}),
                      "multivariate" =  (lm.fit(cbind(1, x), y)$coefficients)[-1]
                      )
    weights <- 1/abs(weights)^nu
    list(weights=weights, nu=nu)
}



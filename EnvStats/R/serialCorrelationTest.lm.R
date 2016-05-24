serialCorrelationTest.lm <-
function (x, test = "rank.von.Neumann", alternative = "two.sided", 
    conf.level = 0.95, ...) 
{
    if (!inherits(x, what = "lm")) 
        stop("x must inherit from the class \"lm\"")
    data.name <- deparse(substitute(x))
    x <- residuals(x)
    test <- match.arg(test, c("rank.von.Neumann", "AR1.yw", "AR1.mle"))
    n <- length(x)
    if (test == "AR1.mle") {
        if (any(is.infinite(x)) || any(is.nan(x))) 
            stop("Residuals from lm fit cannot contain infinite (Inf, -Inf) or undefined (NaN) values")
        bad.obs <- sum(is.na(x))
        if ((n - bad.obs) < 3) 
            stop("Residuals from lm fit must have at least 3 non-missing values.")
    }
    else {
        if (!all(is.finite(x))) 
            stop(paste("Rediduals from lm fit cannot contain missing (NA),", 
                "infinite (Inf, -Inf), or undefined (NaN) values."))
        bad.obs <- 0
        if (n < 3) 
            stop("Residuals from lm fit must have at least 3 values.")
    }
    alternative <- match.arg(alternative, c("two.sided", "greater", 
        "less"))
    ret.list <- serialCorrelationTest.default(x = x, test = test, 
        alternative = alternative, conf.level = conf.level)
    ret.list$parent.of.data <- data.name
    ret.list$data.name <- "Residuals"
    ret.list
}

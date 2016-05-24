eqbinom <-
function (x, size = NULL, p = 0.5, method = "mle/mme/mvue", digits = 0) 
{
    if (!is.vector(p, mode = "numeric") || is.factor(p)) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1.")
    method <- match.arg(method)
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Binomial") 
            stop(paste("'eqbinom' estimates quantiles", "for a binomial distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        size <- x$parameters["size"]
        prob <- x$parameters["prob"]
        n <- x$sample.size
        ret.list <- x
    }
    else {
        if (!((is.vector(x, mode = "numeric") && !is.factor(x)) || 
            is.vector(x, mode = "logical"))) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a numeric or logical vector"))
        data.name <- deparse(substitute(x))
        if (is.null(size)) {
            x <- as.numeric(x)
            if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
                is.not.finite.warning(x)
                x <- x[x.ok]
                warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
            }
            size <- length(x)
            if (size == 0) 
                stop("'x' must contain at least one non-missing value")
            if (!all(x == 0 | x == 1)) 
                stop(paste("When 'size' is not supplied and 'x' is numeric,", 
                  "all non-missing values of 'x' must be 0 or 1."))
            x <- sum(x)
        }
        else {
            if (length(x) != 1 || !is.numeric(x) || !is.finite(x) || 
                x != trunc(x) || x < 0) 
                stop("'x' must be a non-negative integer when 'size' is supplied")
            if (length(size) != 1 || !is.numeric(size) || !is.finite(size) || 
                size != trunc(size) || size < x) 
                stop("'size' must be a postive integer at least as large as 'x'")
            bad.obs <- 0
        }
        ret.list <- ebinom(x, size = size, method = method)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        size <- ret.list$parameters["size"]
        prob <- ret.list$parameters["prob"]
    }
    q <- qbinom(p, size = size, prob = prob)
    if (length(p) == 1 && p == 0.5) 
        names(q) <- "Median"
    else {
        pct <- round(100 * p, digits)
        names(q) <- paste(pct, number.suffix(pct), " %ile", sep = "")
    }
    ret.list <- c(ret.list, list(quantiles = q))
    ret.list$quantile.method <- paste("Quantile(s) Based on\n", 
        space(33), ret.list$method, " Estimators", sep = "")
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}

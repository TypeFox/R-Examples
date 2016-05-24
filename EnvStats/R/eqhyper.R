eqhyper <-
function (x, m = NULL, total = NULL, k = NULL, p = 0.5, method = "mle", 
    digits = 0) 
{
    if (!is.vector(p, mode = "numeric")) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1.")
    method <- match.arg(method, c("mle", "mvue"))
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Hypergeometric") 
            stop(paste("'eqhyper' estimates quantiles", "for a hypergeometric distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        m <- x$parameters["m"]
        n <- x$parameters["n"]
        k <- x$parameters["k"]
        ret.list <- x
    }
    else {
        if (!is.vector(x, mode = "numeric") || length(x) != 1 || 
            !is.finite(x) || x != trunc(x)) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a non-missing integer"))
        bad.obs <- 0
        if (!is.vector(k, mode = "numeric") || length(k) != 1 || 
            !is.finite(k) || k != trunc(k)) 
            stop("Missing values not allowed for 'k' and it must be an integer")
        if ((is.null(m) & is.null(total)) || (!is.null(m) & !is.null(total))) 
            stop("You must supply 'm' or 'total' but not both")
        if (!is.null(m)) {
            if (!is.vector(m, mode = "numeric") || length(m) != 
                1 || !is.finite(m) || m != trunc(m)) 
                stop("Missing values not allowed for 'm' and it must be an integer")
        }
        else {
            if (!is.vector(total, mode = "numeric") || length(total) != 
                1 || !is.finite(total) || total != trunc(total)) 
                stop("Missing values not allowed for 'total' and it must be an integer")
        }
        data.name <- deparse(substitute(x))
        ret.list <- ehyper(x, m = m, total = total, k = k, method = method)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        m <- ret.list$parameters["m"]
        n <- ret.list$parameters["n"]
        k <- ret.list$parameters["k"]
    }
    q <- qhyper(p, m = m, n = n, k = k)
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

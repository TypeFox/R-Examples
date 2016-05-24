ppccNormMultiplyCensoredGofTest <-
function (x, censored, censoring.side = c("left", "right"), distribution = c("norm", 
    "lnorm", "lnormAlt"), est.arg.list = NULL, prob.method = c("hirsch-stedinger", 
    "michael-schucany", "kaplan-meier", "nelson"), plot.pos.con = 0.375) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!((is.vector(censored, mode = "numeric") && !is.factor(censored)) || 
        is.vector(censored, mode = "logical"))) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    data.name <- deparse(substitute(x))
    censoring.name <- deparse(substitute(censored))
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
        is.not.finite.warning(x)
        is.not.finite.warning(as.numeric(censored))
        x <- x[ok]
        censored <- censored[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'censored' removed."))
    }
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    n.cen <- sum(censored)
    if (n.cen == 0) {
        warning(paste("No censored values indicated by 'censored',", 
            "so the function 'gofTest' was called."))
        ret.list <- gofTest(y = x, test = "ppcc", distribution = distribution, 
            est.arg.list = est.arg.list)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        return(ret.list)
    }
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values.")
    censoring.side <- match.arg(censoring.side)
    distribution <- match.arg(distribution)
    if (any(distribution == c("lnorm", "lnormAlt")) && any(x <= 
        0)) 
        stop("All values of 'x' must be positive for a lognormal distribution")
    prob.method <- match.arg(prob.method)
    if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
        1 || plot.pos.con < 0 || plot.pos.con > 1) 
        stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
    x.cen <- x[censored]
    cen.levels <- sort(unique(x.cen))
    K <- length(cen.levels)
    if (K == 1) {
        warning(paste("Only one censoring level indicated by 'censored',", 
            "so the function 'ppccNormSinglyCensoredGofTest' was called."))
        ret.list <- ppccNormSinglyCensoredGofTest(x = x, censoring.side = censoring.side, 
            distribution = distribution, est.arg.list = est.arg.list)
        ret.list$data.name <- data.name
        ret.list$censoring.name <- censoring.name
        ret.list$bad.obs <- bad.obs
        return(ret.list)
    }
    N <- length(x)
    if (N < 5) 
        stop("'x' should contain at least 5 non-missing values.")
    if (N < 20 || N > 5000) 
        warning(paste("Too few or too many observations.  This approximation only works", 
            "if the number of observations is between 20 and 5000."))
    gof.list <- sfMultiplyCensoredGofTest(x = x, censored = censored, 
        censoring.side = censoring.side, distribution = distribution, 
        est.arg.list = est.arg.list, prob.method = prob.method, 
        plot.pos.con = plot.pos.con)
    gof.list$method <- paste("PPCC GOF", "(Multiply Censored Data)", 
        sep = paste("\n", space(33), sep = ""))
    gof.list$data.name <- data.name
    gof.list$censoring.name <- censoring.name
    r <- sqrt(gof.list$statistic)
    names(r) <- "r"
    gof.list$statistic <- r
    gof.list
}

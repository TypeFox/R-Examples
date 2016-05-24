ppccNormSinglyCensoredGofTest <-
function (x, censored, censoring.side = c("left", "right"), distribution = c("norm", 
    "lnorm", "lnormAlt"), est.arg.list = NULL) 
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
    N <- length(x)
    if (N < 5) 
        stop("'x' must contain at least 5 non-missing values.")
    if (N < 20 || N > 5000) 
        warning(paste("Too few or too many observations.  This approximation only works", 
            "if the number of observations is between 20 and 5000."))
    censoring.side <- match.arg(censoring.side)
    distribution <- match.arg(distribution)
    if (any(distribution == c("lnorm", "lnormAlt")) && any(x <= 
        0)) 
        stop("All values of 'x' must be positive for a lognormal distribution")
    T1 <- unique(x[censored])
    if (length(T1) > 1) 
        stop("More than one censoring level.")
    if (censoring.side == "left") {
        if (T1 > min(x.no.cen)) 
            stop(paste("For singly left-censored data,", "all uncensored observations must be bigger than", 
                "or equal to the censoring level."))
    }
    else {
        if (T1 < max(x.no.cen)) 
            stop(paste("For singly right-censored data,", "all uncensored observations must be less than", 
                "or equal to the censoring level."))
    }
    gof.list <- sfSinglyCensoredGofTest(x = x, censored = censored, 
        censoring.side = censoring.side, distribution = distribution, 
        est.arg.list = est.arg.list)
    gof.list$method <- paste("PPCC GOF", "(Singly Censored Data)", 
        sep = paste("\n", space(33), sep = ""))
    gof.list$data.name <- data.name
    gof.list$censoring.name <- censoring.name
    r <- sqrt(gof.list$statistic)
    names(r) <- "r"
    gof.list$statistic <- r
    gof.list
}

gofTestCensored <-
function (x, censored, censoring.side = "left", test = "sf", 
    distribution = "norm", est.arg.list = NULL, prob.method = "hirsch-stedinger", 
    plot.pos.con = 0.375) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if (!is.vector(censored, mode = "numeric") & !is.vector(censored, 
        mode = "logical")) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    if (is.numeric(censored)) {
        index <- is.finite(censored)
        if (!all(is.element(censored[index], 0:1))) 
            stop(paste("When 'censored' is a numeric vector, all non-missing values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
    }
    censoring.name <- deparse(substitute(censored))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    test <- match.arg(test, c("sw", "sf", "ppcc"))
    if (test == "ppcc") 
        test <- "ppccNorm"
    distribution <- match.arg(distribution, c("norm", "lnorm", 
        "lnormAlt"))
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
        is.not.finite.warning(x)
        is.not.finite.warning(as.numeric(censored))
        x <- x[ok]
        censored <- censored[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'censored' removed."))
    }
    if (is.numeric(censored)) 
        censored <- as.logical(censored)
    n.cen <- sum(censored)
    if (n.cen == 0) {
        warning(paste("No censored values indicated by 'censored',", 
            "so the function 'gofTest' was called."))
        ret.list <- gofTest(x = x, test = test, distribution = distribution, 
            est.arg.list = est.arg.list)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        return(ret.list)
    }
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values.")
    if (any(distribution == c("lnorm", "lnormAlt")) && any(x <= 
        0)) 
        stop("All non-missing values of 'x' must be positive for a lognormal distribution")
    multiple <- TRUE
    T.vec <- unique(x[censored])
    if (length(T.vec) == 1) {
        if (censoring.side == "left") {
            if (T.vec <= min(x.no.cen)) 
                multiple <- FALSE
        }
        else {
            if (T.vec >= max(x.no.cen)) 
                multiple <- FALSE
        }
    }
    if (multiple) {
        if (test == "sw") 
            stop(paste("Shapiro-Wilk test not available for multiply censored data.", 
                "Set test='sf' or test='ppcc'."))
        prob.method <- match.arg(prob.method, c("hirsch-stedinger", 
            "michael-schucany", "modified kaplan-meier", "nelson"))
        if (censoring.side == "left" & prob.method == "nelson") 
            stop("Nelson Method not available when censoring.side='left'")
        if (censoring.side == "right" & prob.method == "modified kaplan-meier") 
            stop("Modified Kaplan-Meier Method not available when censoring.side='right'")
        if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
            1 || plot.pos.con < 0 || plot.pos.con > 1) 
            stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
        censoring.type <- "MultiplyCensored"
    }
    else {
        censoring.type <- "SinglyCensored"
    }
    test.name <- paste(test, censoring.type, "GofTest", sep = "")
    arg.list <- switch(test.name, swSinglyCensoredGofTest = , 
        sfSinglyCensoredGofTest = , ppccNormSinglyCensoredGofTest = list(x = x, 
            censored = censored, censoring.side = censoring.side, 
            distribution = distribution, est.arg.list = est.arg.list), 
        sfMultiplyCensoredGofTest = , ppccNormMultiplyCensoredGofTest = list(x = x, 
            censored = censored, censoring.side = censoring.side, 
            distribution = distribution, est.arg.list = est.arg.list, 
            prob.method = prob.method, plot.pos.con = plot.pos.con))
    ret.list <- do.call(test.name, args = arg.list)
    ret.list$data.name <- data.name
    ret.list$censoring.name <- censoring.name
    ret.list
}

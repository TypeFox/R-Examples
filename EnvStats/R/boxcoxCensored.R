boxcoxCensored <-
function (x, censored, censoring.side = "left", lambda = {
    if (optimize) 
        c(-2, 2)
    else seq(-2, 2, by = 0.5)
}, optimize = FALSE, objective.name = "PPCC", eps = .Machine$double.eps, 
    include.x.and.censored = TRUE, prob.method = "michael-schucany", 
    plot.pos.con = 0.375) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!is.vector(censored, mode = "numeric") && !is.vector(censored, 
        mode = "logical")) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    if (!is.vector(lambda, mode = "numeric") || any(!is.finite(lambda))) 
        stop("'lambda' must be a numeric vector with no missing or infinite values")
    if (optimize && length(unique(lambda)) != 2) 
        stop(paste("When optimize=TRUE, 'lambda' must be a vector", 
            "with two unique values that specify the lower and", 
            "upper bounds for the optimization"))
    if (optimize && (1 < min(lambda) || 1 > max(lambda))) 
        stop("When optimize=TRUE, the range of 'lambda' must contain 1")
    data.name <- deparse(substitute(x))
    censoring.name <- deparse(substitute(censored))
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
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
    if (n.cen == 0) 
        stop("No censored values indicated by 'censored'.")
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values.")
    if (any(x <= 0)) 
        stop("All non-missing, finite values of 'x' must be positive")
    objective.name <- match.arg(objective.name, c("PPCC", "Log-Likelihood"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
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
    args.list <- list(x = x, censored = censored, censoring.side = censoring.side, 
        lambda = lambda, optimize = optimize, objective.name = objective.name, 
        eps = eps, include.x.and.censored = include.x.and.censored)
    if (multiple) {
        if (objective.name == "PPCC") {
            prob.method <- match.arg(prob.method, c("michael-schucany", 
                "hirsch-stedinger", "kaplan-meier", "modified kaplan-meier", 
                "nelson"))
            if (!is.vector(plot.pos.con, mode = "numeric") || 
                length(plot.pos.con) != 1 || plot.pos.con < 0 || 
                plot.pos.con > 1) 
                stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
            args.list <- c(args.list, list(prob.method = prob.method, 
                plot.pos.con = plot.pos.con))
        }
        ret.list <- do.call("boxcoxMultiplyCensored", args = args.list)
    }
    else {
        ret.list <- do.call("boxcoxSinglyCensored", args = args.list)
    }
    ret.list$data.name <- data.name
    ret.list$censoring.name <- censoring.name
    ret.list
}

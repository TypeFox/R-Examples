simulateVector <-
function (n, distribution = "norm", param.list = list(mean = 0, 
    sd = 1), sample.method = "SRS", seed = NULL, sorted = FALSE, 
    left.tail.cutoff = ifelse(is.finite(supp.min), 0, .Machine$double.eps), 
    right.tail.cutoff = ifelse(is.finite(supp.max), 0, .Machine$double.eps)) 
{
    if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
        n != trunc(n) || n < 1) 
        stop("'n' must be a positive integer.")
    if (!is.vector(distribution, mode = "character") || length(distribution) != 
        1) 
        stop("'distribution' must be a character string.")
    sample.method <- match.arg(sample.method, c("SRS", "LHS"))
    idist <- charmatch(distribution, c("emp", .Distribution.abb), 
        nomatch = 0)
    if (any(idist == 0)) 
        stop(paste("The argument 'distribution' contains an unknown or", 
            "ambiguous distribution abbreviation. ", "See the help file for 'EnvStats::Distribution.df' for", 
            "a list of distribution abbreviations. ", "You may also use the abbreviation 'emp' to denote", 
            "an empirical distribution"))
    emp <- idist == 1
    if (!is.list(param.list)) 
        stop("'param.list' must be a list.")
    if (!emp) {
        check.da.list <- check.distribution.args.simulate(distribution, 
            param.list, sample.method)
        dist.abb <- check.da.list$dist.abb
        dist.name <- check.da.list$dist.name
        param.list <- check.da.list$param.list
        supp.min <- eval(parse(text = EnvStats::Distribution.df[dist.abb, 
            "Support.Min"]), envir = param.list)
        supp.max <- eval(parse(text = EnvStats::Distribution.df[dist.abb, 
            "Support.Max"]), envir = param.list)
    }
    else {
        dist.abb <- "emp"
        names.param.list <- names(param.list)
        if (length(param.list) < 1 || (sample.method == "SRS" && 
            length(param.list) > 1) || is.null(names.param.list) || 
            any(names.param.list == "") || all(is.na(charmatch(names.param.list, 
            "obs")))) 
            stop(paste("You have specified sampling from ", "an empirical distribution.  ", 
                "For empirical distributions, ", "'param.list' must be a list of the ", 
                "form\n\n\t", "list(obs = ?)\n\n\t", "where ? denotes a numeric vector of empirical ", 
                "observations.  If sample.method=='LHS', this list ", 
                "may also include other arguments to the function ", 
                "qemp(), such as \n\n\t", "list(obs = ?, discrete = TRUE)\n\n\t", 
                "where again ? denotes a numeric vector", sep = ""))
        obs <- param.list$obs
        if ((bad.obs <- sum(!(obs.ok <- is.finite(obs)))) > 0) {
            is.not.finite.warning(obs)
            param.list$obs <- obs[obs.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'obs' removed."))
        }
        if (!is.vector(param.list$obs, mode = "numeric") || length(param.list$obs) <= 
            1) 
            stop(paste("You have specified sampling from ", "an empirical distribution.  ", 
                "In this case, the component of 'param.list' ", 
                "named 'obs' must be a numeric vector ", "with at least two non-missing elements.", 
                sep = ""))
        if (sample.method == "LHS" && length(param.list) > 1 && 
            any(is.na(charmatch(names.param.list, c("obs", "discrete", 
                "prob.method", "plot.pos.con"))))) 
            stop(paste("A component of 'param.list' has a name that is", 
                "an unknown argument to 'qemp'."))
        supp.min <- min(param.list$obs)
        supp.max <- max(param.list$obs)
    }
    if (!is.null(seed)) {
        if (!is.vector(seed, mode = "numeric") || length(seed) != 
            1 || seed != trunc(seed) || seed < 0 || seed > 1000) 
            stop("'seed' must be  an integer between 0 and 1000")
        set.seed(seed)
    }
    if (sample.method == "LHS") {
        if (n < 2) 
            stop(paste("The value of 'n' must be greater than 1", 
                "when sample.method='LHS'"))
        if (!is.vector(left.tail.cutoff, mode = "numeric") || 
            length(left.tail.cutoff) != 1 || left.tail.cutoff < 
            0 || !is.vector(right.tail.cutoff, mode = "numeric") || 
            length(right.tail.cutoff) != 1 || right.tail.cutoff < 
            0 || left.tail.cutoff >= 1 - right.tail.cutoff) 
            stop(paste("'left.tail.cutoff' and 'right.tail.cutoff' must be", 
                "numeric scalars between 0 and 1, and 'left.tail.cutoff'", 
                "must be smaller than 1-'right.tail.cutoff'"))
        if (left.tail.cutoff == 0 && supp.min == -Inf) 
            stop(paste("The value of 'left.tail.cutoff' must be greater", 
                "than 0 for the", dist.name, "distribution since the support", 
                "on the left-hand tail is infinite."))
        if (right.tail.cutoff == 0 && supp.max == Inf) 
            stop(paste("The value of 'right.tail.cutoff' must be greater", 
                "than 0 for the", dist.name, "distribution since the support", 
                "on the right-hand tail is infinite."))
    }
    if (sample.method == "SRS") {
        rname <- paste("r", dist.abb, sep = "")
        x <- do.call(rname, c(list(n = n), param.list))
        if (sorted) 
            x <- sort(x)
    }
    else {
        p <- seq(left.tail.cutoff, 1 - right.tail.cutoff, length = n + 
            1)
        p <- runif(n, min = p[1:n], max = p[2:(n + 1)])
        qname <- paste("q", dist.abb, sep = "")
        x <- do.call(qname, c(list(p = p), param.list))
        if (!sorted) 
            x <- sample(x)
    }
    x
}

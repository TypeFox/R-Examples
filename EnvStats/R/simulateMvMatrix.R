simulateMvMatrix <-
function (n, distributions = c(Var.1 = "norm", Var.2 = "norm"), 
    param.list = list(Var.1 = list(mean = 0, sd = 1), Var.2 = list(mean = 0, 
        sd = 1)), cor.mat = diag(length(distributions)), sample.method = "SRS", 
    seed = NULL, left.tail.cutoff = ifelse(is.finite(supp.min), 
        0, .Machine$double.eps), right.tail.cutoff = ifelse(is.finite(supp.max), 
        0, .Machine$double.eps), tol.1 = .Machine$double.eps, 
    tol.symmetry = .Machine$double.eps, tol.recip.cond.num = .Machine$double.eps, 
    max.iter = 10) 
{
    if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
        n != trunc(n) || n < 1) 
        stop("'n' must be a positive integer.")
    if (!is.vector(distributions, mode = "character")) 
        stop("'distributions' must be a character vector.")
    k <- length(distributions)
    if (k < 1) 
        stop("'distributions' must contain at least one element")
    var.names <- names(distributions)
    lsm <- length(sample.method)
    if (!is.vector(sample.method, mode = "character") || !(lsm == 
        1 || lsm == k)) 
        stop(paste("'sample.method' must be a character vector of length 1", 
            "or of length equal to the number of variables", 
            "(", k, "in this case)."))
    if (any(is.na(charmatch(sample.method, c("LHS", "SRS"))))) 
        stop(paste("All elements of 'sample.method' must be equal to", 
            "'LHS' or 'SRS', or abbreviations of these two character strings"))
    if (lsm == 1) 
        sample.method <- rep(sample.method, k)
    idist <- charmatch(distributions, c("emp", .Distribution.abb), 
        nomatch = 0)
    if (any(idist == 0)) 
        stop(paste("The argument 'distributions' contains unknown or", 
            "ambiguous distribution abbreviations. ", "See the help file for 'EnvStats::Distribution.df' for", 
            "a list of distribution abbreviations. ", "You may also use the abbreviation 'emp' to denote", 
            "an empirical distribution"))
    emp.index <- idist == 1
    if (!is.list(param.list)) 
        stop("'param.list' must be a list.")
    if (length(param.list) != k) 
        stop(paste("The length of 'param.list' must be the same as", 
            "the length of of 'distributions'."))
    names.param.list <- names(param.list)
    if (!is.null(var.names) && !is.null(names.param.list) && 
        any(names.param.list != var.names)) 
        stop(paste("The names associated with 'param.list' must match", 
            "the names associated with 'distributions'."))
    dist.abb <- character(k)
    dist.name <- dist.abb
    supp.min <- numeric(k)
    supp.max <- supp.min
    param.list.tmp <- vector("list", k)
    if (any(!emp.index)) {
        for (i in (1:k)[!emp.index]) {
            check.da.list <- check.distribution.args.simulate(distributions[i], 
                param.list[[i]], sample.method[i])
            dist.abb[i] <- check.da.list$dist.abb
            dist.name[i] <- check.da.list$dist.name
            param.list.tmp[[i]] <- check.da.list$param.list
            supp.min[i] <- eval(parse(text = EnvStats::Distribution.df[dist.abb[i], 
                "Support.Min"]), envir = param.list.tmp[[i]])
            supp.max[i] <- eval(parse(text = EnvStats::Distribution.df[dist.abb[i], 
                "Support.Max"]), envir = param.list.tmp[[i]])
        }
        param.list[!emp.index] <- param.list.tmp[!emp.index]
    }
    if (any(emp.index)) {
        dist.abb[emp.index] <- "emp"
        for (i in (1:k)[emp.index]) {
            param.list.i <- param.list[[i]]
            if (!is.list(param.list.i)) 
                stop("All components of 'param.list' must be lists.")
            names.param.list.i <- names(param.list.i)
            sample.method.i <- sample.method[i]
            if (length(param.list.i) < 1 || (sample.method.i == 
                "SRS" && length(param.list.i) > 1) || is.null(names.param.list.i) || 
                any(names.param.list.i == "") || all(is.na(charmatch(names.param.list.i, 
                "obs")))) 
                stop(paste("You have specified that the ", i, 
                  number.suffix(i), " variable has an empirical distribution.  ", 
                  "For empirical distributions, the corresponding ", 
                  "component of 'param.list' must be a list of the ", 
                  "form\n\n\t", "list(obs = ?)\n\n\t", "where ? denotes a numeric vector of empirical ", 
                  "observations.  If sample.method=='LHS', this list ", 
                  "may also include other arguments to the function ", 
                  "qemp(), such as \n\n\t", "list(obs = ?, discrete = TRUE)\n\n\t", 
                  "where again ? denotes a numeric vector", sep = ""))
            if (!is.vector(param.list.i$obs, mode = "numeric") || 
                length(param.list.i$obs) <= 1) 
                stop(paste("You have specified that the ", i, 
                  number.suffix(i), " variable has an empirical distribution.  ", 
                  "In this case, the corresponding component of ", 
                  "'param.list' must be a list, and within this list ", 
                  "the component called 'obs' must be a numeric vector ", 
                  "with at least two elements.", sep = ""))
            if (sample.method.i == "LHS" && length(param.list.i) > 
                1 && any(is.na(charmatch(names.param.list.i, 
                c("obs", "discrete", "prob.method", "plot.pos.con"))))) 
                stop(paste("A component of the ", i, number.suffix(i), 
                  " component of 'param.list' has a name that is", 
                  "an unknown argument to 'qemp'."))
            supp.min[i] <- min(param.list.i$obs)
            supp.max[i] <- max(param.list.i$obs)
        }
    }
    if (k > 1) {
        if (!is.numeric(cor.mat) || !is.matrix(cor.mat) || !all(dim(cor.mat) == 
            k)) 
            stop(paste("You have specified", k, "distributions, so 'cor.mat'", 
                "must be a", k, "x", k, " numeric matrix"))
        if (any(abs(diag(cor.mat) - 1) > tol.1)) 
            stop(paste("All elements of 'cor.mat' on the diagonal ", 
                "(i.e., all (i,i) elements for i = 1 to", k, 
                ") must be 1.", sep = ""))
        if (!all(cor.mat >= -1 & cor.mat <= 1)) 
            stop("All elements of 'cor.mat' must be between -1 and 1")
        if (any(abs((cor.mat - t(cor.mat))) > tol.symmetry)) 
            stop("'cor.mat' must be symmetric")
        eigen.values <- eigen(cor.mat, symmetric = TRUE, only.values = TRUE)$values
        if (!all(eigen.values > 0) || (min(eigen.values)/max(eigen.values)) < 
            tol.recip.cond.num) 
            stop("'cor.mat' must be a positive definite matrix.")
    }
    if (!is.null(seed)) {
        if (!is.vector(seed, mode = "numeric") || length(seed) != 
            1 || seed != trunc(seed) || seed < 0 || seed > 1000) 
            stop("'seed' must be  an integer between 0 and 1000")
        set.seed(seed)
    }
    if (any(sample.method == "LHS")) {
        if (!is.vector(left.tail.cutoff, mode = "numeric") || 
            length(left.tail.cutoff) != k || any(left.tail.cutoff < 
            0) || !is.vector(right.tail.cutoff, mode = "numeric") || 
            length(right.tail.cutoff) != k || any(right.tail.cutoff < 
            0) || any(left.tail.cutoff >= 1 - right.tail.cutoff)) 
            stop(paste("All values of", "'left.tail.cutoff' and 'right.tail.cutoff' must be", 
                "numeric scalars between 0 and 1, and all values of", 
                "'left.tail.cutoff' must be smaller than the", 
                "corresponding values of 1-'right.tail.cutoff'"))
        if (any(bad.tail.cutoff <- left.tail.cutoff == 0 & supp.min == 
            -Inf)) 
            stop(paste("The value of 'left.tail.cutoff' must be greater", 
                "than 0 for the", (1:k)[bad.tail.cutoff], "distributions since the support", 
                "on the left-hand tail is infinite."))
        if (any(bad.tail.cutoff <- right.tail.cutoff == 0 & supp.max == 
            Inf)) 
            stop(paste("The value of 'right.tail.cutoff' must be greater", 
                "than 0 for the", (1:k)[bad.tail.cutoff], "distributions since the support", 
                "on the right-hand tail is infinite."))
    }
    P <- t(chol(cor.mat))
    R <- matrix(qnorm((1:n)/(n + 1)), n, k)
    cor.R.pos.def <- FALSE
    iter <- 1
    while (!cor.R.pos.def && iter <= max.iter) {
        R <- apply(R, 2, sample)
        cor.R <- cor(R)
        eigen.values <- eigen(cor.R, symmetric = TRUE, only.values = TRUE)$values
        if (all(eigen.values > 0) && (min(eigen.values)/max(eigen.values)) >= 
            tol.recip.cond.num) 
            cor.R.pos.def <- TRUE
        iter <- iter + 1
    }
    if (iter > max.iter) 
        stop(paste("Unable to create an R matrix with a", "positive definite correlation matrix. ", 
            "Increase 'n' and/or increase 'max.iter' and/or", 
            "decrease 'tol.recip.cond.num'."))
    Q <- t(chol(cor.R))
    R <- R %*% t(P %*% solve(Q))
    x <- sapply(1:k, function(i, n, dist.abb, param.list, sample.method, 
        left.tail.cutoff, right.tail.cutoff) {
        simulateVector(n = n, distribution = dist.abb[i], param.list = param.list[[i]], 
            sample.method = sample.method[i], left.tail.cutoff = left.tail.cutoff[i], 
            right.tail.cutoff = right.tail.cutoff[i])
    }, n = n, dist.abb = dist.abb, param.list = param.list, sample.method = sample.method, 
        left.tail.cutoff = left.tail.cutoff, right.tail.cutoff = right.tail.cutoff)
    x <- apply(x, 2, sort)
    x <- sapply(1:k, function(i, x, R) {
        x[rank(R[, i]), i]
    }, x = x, R = R)
    if (!is.null(var.names)) 
        dimnames(x) <- list(NULL, var.names)
    x
}

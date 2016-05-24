ebinom <-
function (x, size = NULL, method = "mle/mme/mvue", ci = FALSE, 
    ci.type = "two-sided", ci.method = "score", correct = TRUE, 
    var.denom = "n", conf.level = 0.95, warn = TRUE) 
{
    if (x.fac <- is.factor(x)) {
        if (length(levels(x)) != 2) 
            stop("When 'x' is a factor it must have 2 levels")
    }
    else {
        if (!is.vector(x, mode = "numeric") && !is.vector(x, 
            mode = "logical")) 
            stop("'x' must be a numeric or logical vector or else a factor with two levels")
    }
    data.name <- deparse(substitute(x))
    method <- match.arg(method)
    if (is.null(size)) {
        x <- as.numeric(x)
        if (x.fac) 
            x <- x - 1
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            x <- x[x.ok]
            string <- ifelse(bad.obs > 1, "observations", "observation")
            warning(paste(bad.obs, string, "with NA/NaN/Inf in 'x' removed."))
        }
        size <- length(x)
        if (size == 0) 
            stop("'x' must contain at least one non-missing value")
        if (!all(x == 0 | x == 1)) 
            stop(paste("When 'size' is not supplied and 'x' is numeric,", 
                "all non-missing values of 'x' must be 0 or 1."))
        prob <- mean(x)
    }
    else {
        if (length(x) != 1 || !is.numeric(x) || !is.finite(x) || 
            x != trunc(x) || x < 0) 
            stop("'x' must be a non-negative integer when 'size' is supplied")
        if (length(size) != 1 || !is.numeric(size) || !is.finite(size) || 
            size != trunc(size) || size < x) 
            stop("'size' must be a postive integer at least as large as 'x'")
        bad.obs <- 0
        prob <- x/size
    }
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method, c("score", "exact", 
            "adjusted Wald", "Wald"))
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        switch(ci.method, score = {
            ci.obj <- ci.binom.score(x = prob * size, size = size, 
                alpha = (1 - conf.level), ci.type = ci.type, 
                correct = correct)
        }, exact = {
            ci.obj <- ci.binom.exact(x = prob * size, size = size, 
                alpha = (1 - conf.level), ci.type = ci.type)
        }, `adjusted Wald` = {
            if (prob == 0 || prob == 1) stop(paste("All successes or all failures. ", 
                "You must use the 'score' or 'exact' method for this case."))
            ci.obj <- ci.binom.adjusted.Wald(x = prob * size, 
                size = size, alpha = (1 - conf.level), ci.type = ci.type)
        }, Wald = {
            if (prob == 0 || prob == 1) stop(paste("All successes or all failures. ", 
                "You must use the 'score' or 'exact' method for this case."))
            var.denom <- match.arg(var.denom, c("n", "n-1"))
            if (var.denom == "n-1" && size < 2) stop("'size' must be greater than 1 to use the 'n-1' variance estimator")
            ci.obj <- ci.binom.Wald(x = prob * size, size = size, 
                alpha = (1 - conf.level), ci.type = ci.type, 
                correct = correct, var.denom = var.denom, warn = warn)
        })
    }
    dist.params <- c(size = size, prob = prob)
    ret.list <- list(distribution = "Binomial", sample.size = size, 
        parameters = dist.params, n.param.est = 1, method = paste(method, 
            "for 'prob'"), data.name = data.name, bad.obs = bad.obs)
    if (ci) 
        ret.list <- c(ret.list, list(interval = ci.obj))
    oldClass(ret.list) <- "estimate"
    ret.list
}

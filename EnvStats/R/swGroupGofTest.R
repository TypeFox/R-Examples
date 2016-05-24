swGroupGofTest <-
function (x, group, distribution = c("norm", "lnorm", "lnormAlt", 
    "lnorm3", "zmnorm", "zmlnorm", "zmlnormAlt")) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    grouping.variable <- deparse(substitute(group))
    group <- unclass(group)
    if (any(!is.finite(group))) 
        stop("NA's/Inf's not allowed in 'group'.")
    group <- factor(group)
    if (length(x) != length(group)) 
        stop("'x' and 'group' must have the same length.")
    x.list <- split(x, group)
    n.grps <- length(x.list)
    c.names <- names(x.list)
    if (is.null(c.names)) {
        c.names <- paste("V", 1:n.grps, sep = ".")
    }
    if (n.grps == 1) {
        ret.list <- swGofTest(unlist(x.list), distribution = distribution)
        ret.list$data.name <- data.name[1]
        warning("Only one group supplied, so the function 'swGofTest' was called.")
        return(ret.list)
    }
    distribution <- match.arg(distribution)
    test.fcn <- function(x, distribution) {
        x.ok <- x[is.finite(x)]
        if (any(distribution == c("lnorm", "lnormAlt")) && any(x.ok <= 
            0)) 
            stop("All values must be positive for a lognormal distribution")
        if (any(distribution == c("zmlnorm", "zmlnormAlt")) && 
            any(x.ok < 0)) 
            stop(paste("All values must be non-negative for a", 
                "zero-modified lognormal distribution"))
        n <- switch(distribution, zmnorm = sum(x.ok != 0), zmlnorm = , 
            zmlnormAlt = sum(x.ok > 0), length(x.ok))
        if (distribution == "lnorm3") {
            if (n < 5 || length(unique(x.ok)) < 3) 
                stop(paste("Each group of observations must contain", 
                  "at least 5 non-missing values,", "and at least 3 distinct values."))
            if (n > 2000) 
                stop(paste("Too many observations in a group. ", 
                  "This approximation only works", "if the number of observations within any group", 
                  "is between 5 and 2000"))
        }
        else {
            if (n < 4 || all(x.ok == x.ok[1])) 
                stop(paste("Each group of observations must contain", 
                  "at least 4 non-missing values,", "and at least 2 distinct values."))
            if (n > 5000) 
                stop(paste("Too many observations in a group. ", 
                  "This approximation only works", "if the number of observations is between 4 and 5000"))
        }
    }
    sapply(x.list, test.fcn, distribution = distribution)
    sw.fcn <- function(x, distribution) {
        swGofTest(x, distribution = distribution)[c("sample.size", 
            "p.value", "bad.obs")]
    }
    dum.list <- lapply(x.list, sw.fcn, distribution = distribution)
    sample.size <- sapply(dum.list, function(x) x$sample.size)
    names(sample.size) <- c.names
    p.value.vec <- sapply(dum.list, function(x) x$p.value)
    names(p.value.vec) <- c.names
    bad.obs.vec <- sapply(dum.list, function(x) x$bad.obs)
    names(bad.obs.vec) <- c.names
    G.vec <- qnorm(p.value.vec)
    G <- sum(G.vec)/sqrt(n.grps)
    C.vec <- -2 * log(p.value.vec)
    C <- sum(C.vec)
    df <- 2 * n.grps
    p.values <- c(pnorm(G), 1 - pchisq(C, df = df))
    names(p.values) <- c("z (G)", "Chi-Square (C)")
    dist.abb <- distribution
    distribution <- .Distribution.name[match(dist.abb, .Distribution.abb)]
    ret.list <- list(distribution = distribution, dist.abb = dist.abb, 
        statistic = c(`z (G)` = G, `Chi-Square (C)` = C), sample.size = sample.size, 
        parameters = c(df = df), p.value = c(p.value.vec, p.values), 
        alternative = paste("At least one group", "does not come from a", 
            paste(distribution, "Distribution."), sep = paste("\n", 
                space(33), sep = "")), method = "Shapiro-Wilk Group Test", 
        data.name = data.name, grouping.variable = grouping.variable, 
        bad.obs = bad.obs.vec, n.groups = n.grps, group.names = c.names, 
        G.vec = G.vec, C.vec = C.vec)
    oldClass(ret.list) <- "gofGroup"
    ret.list
}

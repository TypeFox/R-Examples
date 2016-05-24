boxplotlog <-
function (x, ..., range = 1.5, width = NULL, varwidth = FALSE,
    notch = FALSE, outline = TRUE, names, plot = TRUE, border = par("fg"),
    col = NULL, log = "", pars = list(boxwex = 0.8, staplewex = 0.5,
        outwex = 0.5), horizontal = FALSE, add = FALSE, at = NULL)
{
# draw log-boxplot:

    args <- list(x, ...)
    namedargs <- if (!is.null(attributes(args)$names))
        attributes(args)$names != ""
    else rep(FALSE, length.out = length(args))
    pars <- c(args[namedargs], pars)
    groups <- if (is.list(x))
        x
    else args[!namedargs]
if (is.list(x))
    groups.l <- log10(x)
else {
    groups.l <- args[!namedargs]
    #groups.l[[1]] <- log10(groups.l[[1]])
      for (i in 1:length(groups.l)){
        groups.l[[i]] <- log10(groups.l[[i]])
      }
}

    if (0 == (n <- length(groups)))
        stop("invalid first argument")
    if (length(class(groups)))
        groups <- unclass(groups)
    if (!missing(names))
        attr(groups, "names") <- names
    else {
        if (is.null(attr(groups, "names")))
            attr(groups, "names") <- 1:n
        names <- attr(groups, "names")
    }
    if (0 == (n <- length(groups.l)))
        stop("invalid first argument")
    if (length(class(groups.l)))
        groups.l <- unclass(groups.l)
    if (!missing(names))
        attr(groups.l, "names") <- names
    else {
        if (is.null(attr(groups.l, "names")))
            attr(groups.l, "names") <- 1:n
        names <- attr(groups.l, "names")
    }
    for (i in 1:n) groups.l[i] <- list(boxplot.stats(groups.l[[i]], range))
    stats.l <- matrix(0, nrow = 5, ncol = n)
    conf.l <- matrix(0, nrow = 2, ncol = n)
    ng.l <- out.l <- group.l <- numeric(0)
    ct.l <- 1
    for (i in groups.l) {
        stats.l[, ct.l] <- i$stats
        conf.l[, ct.l] <- i$conf
        ng.l <- c(ng.l, i$n)
        if ((lo.l <- length(i$out))) {
           out.l <- c(out.l, i$out)
            group.l <- c(group.l, rep.int(ct.l, lo.l))
        }
        ct.l <- ct.l + 1
    }
    for (i in 1:n) groups[i] <- list(boxplot.stats(groups[[i]], range))
    stats <- matrix(0, nrow = 5, ncol = n)
    conf <- matrix(0, nrow = 2, ncol = n)
    stats <- matrix(0, nrow = 5, ncol = n)
    conf <- matrix(0, nrow = 2, ncol = n)
    ng <- out <- group <- numeric(0)
    ct <- 1
    for (i in groups) {
        stats[, ct] <- i$stats
        stats[c(1,5), ct] <- 10^stats.l[c(1,5), ct]
        conf[, ct] <- i$conf
        ng <- c(ng, i$n)
        if ((lo <- length(i$out))) {
            out <- c(out, i$out)
            group <- c(group, rep.int(ct, lo))
        }
        ct <- ct + 1
    }
    z <- list(stats = stats, n = ng, conf = conf, out = out,
        group = group, names = names)
    z.l <- list(stats = stats, n = ng, conf = conf, out = 10^out.l,
        group = group.l, names = names)
    if (plot) {
        bxp(z.l, width, varwidth = varwidth, notch = notch, log = log,
            border = border, col = col, pars = pars, outline = outline,
            horizontal = horizontal, add = add, at = at)
        invisible(z.l)
    }
    else z.l
}

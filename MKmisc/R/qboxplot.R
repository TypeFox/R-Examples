## Boxplots based on function quantile
qboxplot <- function(x, ...){
  UseMethod('qboxplot')
}

qbxp.stats <- function(x, coef = 1.5, do.conf = TRUE, do.out = TRUE, 
                           type = 7){
    if (coef < 0) stop("'coef' must not be negative")
    nna <- !is.na(x)
    n <- sum(nna)
    stats <- fiveNS(x, na.rm = TRUE, type = type)
    iqr <- diff(stats[c(2, 4)])
    if (coef == 0) 
        do.out <- FALSE
    else {
        out <- if (!is.na(iqr)) {
            x < (stats[2] - coef * iqr) | x > (stats[4] + coef * 
                iqr)
        }
        else !is.finite(x)
        if (any(out[nna], na.rm = TRUE)) 
            stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
    }
    conf <- if (do.conf) 
        stats[3] + c(-1.58, 1.58) * iqr/sqrt(n)
    list(stats = stats, n = n, conf = conf, out = if (do.out) x[out & 
        nna] else numeric(0))
}

qboxplot.default <- function(x, ..., range = 1.5, width = NULL, varwidth = FALSE, 
    notch = FALSE, outline = TRUE, names, plot = TRUE, border = par("fg"), 
    col = NULL, log = "", pars = list(boxwex = 0.8, staplewex = 0.5, 
        outwex = 0.5), horizontal = FALSE, add = FALSE, at = NULL, type = 7){
    args <- list(x, ...)
    namedargs <- if (!is.null(attributes(args)$names)) 
        attributes(args)$names != ""
    else rep(FALSE, length.out = length(args))
    groups <- if (is.list(x)) 
        x
    else args[!namedargs]
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
    cls <- sapply(groups, function(x) class(x)[1])
    cl <- if (all(cls == cls[1])) 
        cls[1]
    else NULL
    for (i in 1:n) groups[i] <- list(qbxp.stats(unclass(groups[[i]]), 
        coef = range, type = type))
    stats <- matrix(0, nrow = 5, ncol = n)
    conf <- matrix(0, nrow = 2, ncol = n)
    ng <- out <- group <- numeric(0)
    ct <- 1
    for (i in groups) {
        stats[, ct] <- i$stats
        conf[, ct] <- i$conf
        ng <- c(ng, i$n)
        if ((lo <- length(i$out))) {
            out <- c(out, i$out)
            group <- c(group, rep.int(ct, lo))
        }
        ct <- ct + 1
    }
    if (length(cl) && cl != "numeric") 
        oldClass(stats) <- cl
    z <- list(stats = stats, n = ng, conf = conf, out = out, 
        group = group, names = names)
    if (plot) {
        if (is.null(pars$boxfill) && is.null(args$boxfill)) 
            pars$boxfill <- col
        do.call("bxp", c(list(z, notch = notch, width = width, 
            varwidth = varwidth, log = log, border = border, 
            pars = pars, outline = outline, horizontal = horizontal, 
            add = add, at = at), args[namedargs]))
        invisible(z)
    }
    else z
}

qboxplot.formula <- function(formula, data = NULL, ..., subset, na.action = NULL, type = 7){
    if (missing(formula) || (length(formula) != 3)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action
    m$type <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    qboxplot(split(mf[[response]], mf[-response]), ..., type = type)
}

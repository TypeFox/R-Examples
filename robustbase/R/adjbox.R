#### Skewness (MC) - Adjusted Boxplots

### modeled closely after  boxplot() etc in  R/src/library/graphics/R/boxplot.R :

adjbox <- function(x, ...) UseMethod("adjbox")

adjbox.default <- function (x, ..., range = 1.5, doReflect=FALSE, width = NULL, varwidth = FALSE,
    notch = FALSE, outline = TRUE, names, plot = TRUE, border = par("fg"),
    col = NULL, log = "", pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
			    horizontal = FALSE, add = FALSE, at = NULL)
{
    args <- list(x, ...)
    namedargs <-
	if(!is.null(attributes(args)$names))
	    attributes(args)$names != ""
	else
	    logical(length(args))# all FALSE
    ## pars <- c(args[namedargs], pars)
    groups <- if(is.list(x)) x else args[!namedargs]
    if(0 == (n <- length(groups)))
	stop("invalid first argument")
    if(length(class(groups)))
	groups <- unclass(groups)
    if(!missing(names))
	attr(groups, "names") <- names
    else {
	if(is.null(attr(groups, "names")))
	    attr(groups, "names") <- 1:n
	names <- attr(groups, "names")
    }
    cls <- sapply(groups, function(x) class(x)[1])
    cl <- if(all(cls == cls[1])) cls[1] # else NULL
    for (i in 1:n)
	groups[i] <- list(adjboxStats(unclass(groups[[i]]),
				      coef=range, doReflect=doReflect)) # do.conf=notch)
    stats <- matrix(0, nrow=5, ncol=n)
    conf <- fence <- matrix(0, nrow=2, ncol=n)
    ng <- out <- group <- numeric(0)
    ct <- 1
    for(i in groups) {
	stats[,ct] <- i$stats
	conf [,ct] <- i$conf
	fence[,ct] <- i$fence
	ng <- c(ng, i$n)
	if((lo <- length(i$out))) {
	    out	  <- c(out,i$out)
	    group <- c(group, rep.int(ct, lo))
	}
	ct <- ct+1
    }
    if(length(cl) && cl != "numeric") oldClass(stats) <- cl
    z <- list(stats = stats, n = ng, conf = conf, fence = fence,
              out = out, group = group, names = names)
    if(plot) {
        if(is.null(pars$boxfill) && is.null(args$boxfill)) pars$boxfill <- col
        do.call("bxp",
                c(list(z, notch = notch, width = width, varwidth = varwidth,
                       log = log, border = border, pars = pars,
                       outline = outline, horizontal = horizontal, add = add,
                       at = at), args[namedargs]))
	invisible(z)
    }
    else z
}


adjbox.formula <- function (formula, data = NULL, ..., subset, na.action = NULL)
{
    if(missing(formula) || (length(formula) != 3))
	stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
	m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action # force use of default for this method
    ## require(stats, quietly = TRUE): model.frame
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    adjbox(split(mf[[response]], mf[-response]), ...)
}


## modeled after boxplot.stats()   from	 R/src/library/grDevices/R/calc.R :
adjboxStats <- function(x, coef = 1.5, a = -4, b = 3,
			do.conf = TRUE, do.out = TRUE, ...)
{
    if(coef < 0) stop("'coef' must not be negative")
    nna <- !is.na(x)
    n <- sum(nna)# including +/- Inf
    stats <- fivenum(x, na.rm = TRUE)
    iqr <- diff(stats[c(2, 4)])
    fence <- rep(NA_real_, 2)
    if(coef == 0)
	do.out <- FALSE # no whiskers to be drawn
    else { ## coef > 0
	out <- if (!is.na(iqr)) {
	    medc <- mc(x, ..., na.rm = TRUE)
	    fence <-
		if (medc >= 0)
		    c(stats[2] - coef * exp(a * medc) * iqr,
		      stats[4] + coef * exp(b * medc) * iqr)
		else
		    c(stats[2] - coef * exp(-b * medc) * iqr,
		      stats[4] + coef * exp(-a * medc) * iqr)

	    x < fence[1] | fence[2] < x
	}
	else !is.finite(x)
	if (any(out[nna], na.rm = TRUE))
	    stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
    }
    conf <- if (do.conf) stats[3] + c(-1.58, 1.58) * iqr/sqrt(n)
    list(stats = stats, n = n, conf = conf, fence = fence,
	 out = if (do.out) x[out & nna] else numeric(0))
}

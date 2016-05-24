# checked in 26 December 2009 by J. Fox
# 2012-12-12: Fixed Boxplot.default() so that it works properly when g is numeric. J. Fox
# 2013-04-10: handles at argument properly, contribution of Steve Ellison. J. Fox

# 2013-08-19: removed loading of stats package. J. Fox

Boxplot <- function(y, ...){
	UseMethod("Boxplot")
}

Boxplot.default <- function (y, g, labels, id.method = c("y", "identify", "none"), 
    id.n = 10, xlab, ylab, ...) {
    id.method <- match.arg(id.method)
    if (missing(ylab)) 
        ylab <- deparse(substitute(y))
    if (missing(labels)) 
        labels <- seq(along = y)
    pars <- list(...)
    if (missing(g)) {
        valid <- complete.cases(y, labels)
        y <- y[valid]
        labels <- labels[valid]
        b <- boxplot(y, ylab = ylab, ...)
        if (id.method == "none" | id.n == 0) 
            return(invisible(NULL))
        else if (id.method == "identify") {
            res <- identify(rep(1, length(y)), y, labels)
            return(if (length(res) == 0) invisible(NULL) else labels[res])
        }
        else if (length(b$out) > 0) {
            sel <- y %in% b$out
            yy <- y[sel]
            labs <- labels[sel]
            which.low <- yy < b$stats[1, 1]
            y.low <- yy[which.low]
            labs.low <- labs[which.low]
            if (length(y.low) > id.n) {
                ord.low <- order(y.low)[1:id.n]
                y.low <- y.low[ord.low]
                labs.low <- labs.low[ord.low]
            }
            which.high <- yy > b$stats[5, 1]
            y.high <- yy[which.high]
            labs.high <- labs[which.high]
            if (length(y.high) > id.n) {
                ord.high <- order(y.high, decreasing = TRUE)[1:id.n]
                y.high <- y.high[ord.high]
                labs.high <- labs.high[ord.high]
            }
            labs <- c(labs.low, labs.high)
            at <- if(!is.null(pars$at)) pars$at else 1			#@@@
            text(at, c(y.low, y.high), labs, pos = 2)			#@@@
            return(if (length(labs) == 0) invisible(NULL) else labs)
        }
        else return(invisible(NULL))
    }
    else {
        if (missing(xlab)) 
            xlab = deparse(substitute(g))
        valid <- complete.cases(y, labels, g)
        y <- y[valid]
        labels <- labels[valid]
        g <- g[valid]
        b <- boxplot(split(y, g), ylab = ylab, xlab = xlab, ...)
        levels <- if (is.factor(g)) 
            levels(g)
        else sort(unique(g))
        gg <- as.numeric(g)
        if (id.method == "none" | id.n == 0) 
            return(invisible(NULL))
        else if (id.method == "identify") {
            res <- identify(gg, y, labels)
            return(if (length(res) == 0) invisible(NULL) else labels[res])
        }
        else {
            midx <- mean(par("usr")[1:2])
            identified <- character(0)
            if (length(b$out) > 0) {
                groups <- unique(b$group)
                for (group in groups) {
                    grp <- g == levels[group]
                    yy <- y[grp]
                    labs <- labels[grp]
                    sel <- yy %in% b$out[b$group == group]
                    yy <- yy[sel]
                    glabs <- labs[sel]
                    which.low <- yy < b$stats[1, group]
                    y.low <- yy[which.low]
                    labs.low <- glabs[which.low]
                    if (length(y.low) > id.n) {
                        ord.low <- order(y.low)[1:id.n]
                        y.low <- y.low[ord.low]
                        labs.low <- labs.low[ord.low]
                    }
                    which.high <- yy > b$stats[5, group]
                    y.high <- yy[which.high]
                    labs.high <- glabs[which.high]
                    if (length(y.high) > id.n) {
                        ord.high <- order(y.high, decreasing = TRUE)[1:id.n]
                        y.high <- y.high[ord.high]
                        labs.high <- labs.high[ord.high]
                    }
                    pos <- if (group < midx) 
                        4
                    else 2
                    at <- if(!is.null(pars$at)) pars$at[group] else group
                    text(at, c(y.low, y.high), c(labs.low, labs.high),
                        pos = pos)
                    identified <- c(identified, c(labs.low, labs.high))
                }
            }
            return(if (length(identified) == 0) invisible(NULL) else identified)
        }
    }
}


Boxplot.formula <- function(formula, data=NULL, subset, na.action=NULL, labels., 
  id.method=c("y", "identify", "none"), xlab, ylab, ...){
	# much of this function adapted from graphics:boxplot.formula
	id.method <- match.arg(id.method)
	m <- match.call(expand.dots = FALSE)
	if (is.matrix(eval(m$data, parent.frame()))) 
		m$data <- as.data.frame(data)
	m$xlab <- m$ylab <- m$id.method <- m$... <- NULL
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")
	mf <- eval(m, parent.frame())
	if (missing(labels.)) mf$"(labels.)" <- rownames(mf)
	lab.var <- which(names(mf) == "(labels.)")
	if (length(formula) == 3){
		response <- attr(attr(mf, "terms"), "response")
		if (missing(ylab)) ylab <- names(mf)[response]
		if (missing(xlab)) xlab <- names(mf)[-c(response, lab.var)]
        x <- mf[, -c(response, lab.var)]
		if (is.data.frame(x)) x <- do.call("interaction", as.list(x))
        if (length(xlab) > 1) xlab <- paste(xlab, collapse="*")
		Boxplot(mf[[response]], x, labels=mf[[lab.var]], 
            xlab=xlab, ylab=ylab, id.method=id.method, ...)
	}
	else if (length(formula) == 2){
		if (missing(ylab)) ylab <- names(mf)[-lab.var]
		Boxplot(mf[, -lab.var], labels=mf[[lab.var]], ylab=ylab, id.method=id.method, ...)
	}
	else stop("improper Boxplot formula")   
}

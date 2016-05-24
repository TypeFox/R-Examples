# fancy scatterplots  (J. Fox)

# 2010-09-05: J. Fox: changed color choice
# 2010-09-16: fixed point color when col is length 1
# 2010-12-19: J. Fox: added argument legend.coords to place legend.
# 2011-01-15: J. Fox: If x is a factor, calls Boxplot()
# 2011-03-08: J. Fox: changed col argument
# 2012-04-18: J. Fox: fixed labels argument in scatterplot.formula().
# 2012-04-24: J. Fox: further fix to labels
# 2012-09-12: J. Fox: modified treatment of smoother; added loessLine(), gamLine(), quantregLine(). 
# 2012-09-17: S. Weisberg:  smoothers moved to scatterplotSmoothers.R, defaults changed
# 2012-09-19: J. Fox: restored smooth and span arguments for backwards compatibility
# 2013-02-07: S. Weisberg:  modifed call to showLabels to work correctly with groups
# 2014-09-04: J. Fox: empty groups produce warning rather than error
# 2015-07-17: J. Fox: improved above-plot legends.
# 2015-08-05: J. Fox: fixed sp()

reg <- function(reg.line, x, y, col, lwd, lty, log.x, log.y){
    if(log.x) x <- log(x)
    if(log.y) y <- log(y)
    mod <- reg.line(y ~ x)
    y.hat <- fitted.values(mod)
    x <- model.matrix(mod)[, 2]
    min <- which.min(x)
    max <- which.max(x)
    if (!log.x){
        x1 <- x[min]
        x2 <- x[max]
    }
    else {
        x1 <- exp(x[min])
        x2 <- exp(x[max])
    }
    if (!log.y){
        y1 <- y.hat[min]
        y2 <- y.hat[max]
    }
    else {
        y1 <- exp(y.hat[min])
        y2 <- exp(y.hat[max])
    }
    lines(c(x1, x2), c(y1, y2), lwd=lwd, col=col, lty=lty)
}

find.legend.columns <- function(n, target=min(4, n)){
  rem <- n %% target
  if (rem != 0 && rem < target/2) target <- target - 1
  target
}

scatterplot <- function(x, ...){
    UseMethod("scatterplot", x)
}

scatterplot.formula <- function (formula, data, subset, xlab, ylab, legend.title, legend.coords, labels, ...) {
    na.save <- options(na.action=na.omit)
    on.exit(options(na.save))
    na.pass <- function(dframe) dframe
    m <- match.call(expand.dots=FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
        m$data <- as.data.frame(data)
    m$na.action <- na.pass
    m$legend.coords <- m$legend.title <- m$labels <- m$xlab <- m$ylab <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    if (!inherits(formula, "formula") | length(formula) != 3) 
        stop("invalid formula")    
    formula <- as.character(c(formula))
    formula <- as.formula(sub("\\|", "+", formula))
    m$formula <- formula
    if (missing(data)){ 
        X <- na.omit(eval(m, parent.frame()))
        if (missing(labels)) labels <- gsub("X", "", row.names(X)) 
    }
    else{
        X <- eval(m, parent.frame())
        if (missing(labels)) labels <- row.names(X)
    }
    names <- names(X)
    if (missing(xlab)) xlab <- names[2]
    if (missing(ylab)) ylab <- names[1]
    if (ncol(X) == 2) scatterplot(X[,2], X[,1], xlab=xlab, ylab=ylab, 
                                  labels=labels, ...)
    else {
        if (missing(legend.title)) legend.title <- names[3]
        scatterplot(X[,2], X[,1], groups=X[,3], xlab=xlab, ylab=ylab,  
                    legend.title=legend.title, legend.coords=legend.coords, labels=labels, ...)
    }
}

scatterplot.default <- function(x, y, smoother=loessLine, smoother.args=list(), smooth, span,
                                spread=!by.groups, reg.line=lm,
                                boxplots=if (by.groups) "" else "xy",
                                xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), las=par("las"),
                                lwd=1, lty=1,
                                labels, id.method = "mahal", 
                                id.n = if(id.method[1]=="identify") length(x) else 0, 
                                id.cex = 1, id.col = palette()[1],
                                log="", jitter=list(), xlim=NULL, ylim=NULL,
                                cex=par("cex"), cex.axis=par("cex.axis"), cex.lab=par("cex.lab"), 
                                cex.main=par("cex.main"), cex.sub=par("cex.sub"), 
                                groups, by.groups=!missing(groups), 
                                legend.title=deparse(substitute(groups)), legend.coords, legend.columns,
                                ellipse=FALSE, levels=c(.5, .95), robust=TRUE,
                                col=if (n.groups == 1) palette()[3:1] else rep(palette(), length=n.groups),
                                pch=1:n.groups, 
                                legend.plot=!missing(groups), reset.par=TRUE, grid=TRUE, ...){
    logged <- function(axis=c("x", "y")){
        axis <- match.arg(axis)
        0 != length(grep(axis, log))
    }
    hbox <- function(x){
        if (logged("x")){
            log.x <- "x"
            .x <- log(x)        
        }
        else {
            log.x <- ""
            .x <- x
        }
        plot(x, seq(0, 1, length=length(x)), type="n", axes=FALSE, xlab="", ylab="", log=log.x, xlim=xlim)
        res <- boxplot.stats(.x, coef = 1.5, do.conf=FALSE)
        if (logged("x")){
            res$stats <- exp(res$stats)
            if (!is.null(res$out)) res$out <- exp(res$out)
        }
        LW <- res$stats[1]
        Q1 <- res$stats[2]
        M <- res$stats[3]
        Q3 <- res$stats[4]
        UW <- res$stats[5]
        lines(c(Q1, Q1, Q3, Q3, Q1), c(0, 1, 1, 0, 0))
        lines(c(M, M), c(0, 1))
        lines(c(LW, Q1), c(.5, .5))
        lines(c(Q3, UW), c(.5, .5))
        if (!is.null(res$out)) points(res$out, rep(.5, length(res$out)), cex=cex)
    }
    vbox <- function(y){
        if (logged("y")){
            log.y <- "y"
            .y <- log(y)
        }
        else {
            log.y <- ""
            .y <- y
        }
        plot(seq(0, 1, length=length(y)), y, type="n", axes=FALSE, xlab="", ylab="", log=log.y, ylim=ylim)
        res <- boxplot.stats(.y, coef = 1.5, do.conf=FALSE)
        if (logged("y")){
            res$stats <- exp(res$stats)
            if (!is.null(res$out)) res$out <- exp(res$out)
        }
        LW <- res$stats[1]
        Q1 <- res$stats[2]
        M <- res$stats[3]
        Q3 <- res$stats[4]
        UW <- res$stats[5]
        lines(c(0, 1, 1, 0, 0), c(Q1, Q1, Q3, Q3, Q1))
        lines(c(0, 1), c(M, M))
        lines(c(.5, .5), c(LW, Q1))
        lines(c(.5, .5), c(Q3, UW))
        if (!is.null(res$out)) points(rep(.5, length(res$out)), res$out, cex=cex)
    }
    # force evaluation of some arguments
    by.groups
    legend.plot
    legend.title
    # smooth and span for backwards compatibility
    if (!missing(smooth)) {
        smoother <- if (isTRUE(smooth)) loessLine else FALSE
    }
    if (!missing(span)) smoother.args$span <- span
    if (is.character(family)) family <- eval(parse(text=family))
    if (missing(labels)){
        labels <- if (is.null(names(y)))
            seq(along=y)
        else names(y)
    }
    if (length(labels) != length(y)) stop("labels argument is the wrong length")
    if (is.factor(x)) {
        if (!(id.method %in% c("y", "identify", "none"))) id.method <- "y"
        return(Boxplot(y, x, id.method="y", labels=labels, xlab=xlab, ylab=ylab))
    }
    mar <- par("mar")
    mfcol <- par("mfcol")
    if (reset.par) on.exit(par(mar=mar, mfcol=mfcol))
    if( FALSE == boxplots) boxplots <- ""
    if (!missing(groups)){
        data <- na.omit(data.frame(groups, x, y, labels, stringsAsFactors=FALSE))
        groups <- data[,1]
        if (!is.factor(groups)) groups <- as.factor(groups)
        .x <- data[,2]
        .y <- data[,3]
        labels <- data[,4]
        top <- if (legend.plot && missing(legend.coords)){
            if (missing(legend.columns)) legend.columns <- find.legend.columns(nlevels(groups))
            4 + ceiling(nlevels(groups))/legend.columns
        }
        else mar[3]
        if (legend.plot && !missing(legend.coords) && missing(legend.columns)){
          legend.columns <- 1
        }
    }
    else {
        .x <- x
        .y <- y
        top <- mar[3]
        groups <- factor(rep(1, length(.x)))
    }
    xbox <- length(grep("x", boxplots)) > 0
    ybox <- length(grep("y", boxplots)) > 0
    if (xbox && ybox)
        layout(matrix(c(1, 0, 3, 2), 2, 2),
               widths = c(5, 95),
               heights= c(95, 5))
    else if (ybox)
        layout(matrix(c(1, 2),1, 2),
               widths = c(5, 95),
               heights= 100)
    else if (xbox)
        layout(matrix(c(2, 1), 2, 1),
               widths = 100,
               heights= c(95, 5))
    else layout (matrix(1, 1, 1),
                 widths=100, heights=100)
    par(mar=c(mar[1], 0, top, 0))
    if (ybox > 0) vbox(.y) 
    par(mar=c(0, mar[2], 0, mar[4]))
    if (xbox > 0) hbox(.x) 
    par(mar=c(mar[1:2], top, mar[4]))
    plot(.x, .y, xlab=xlab, ylab=ylab, las=las, log=log, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab,
         cex.main=cex.main, cex.sub=cex.sub, type="n", xlim=xlim, ylim=ylim, ...)
    if(grid){
        grid(lty=1, equilogs=FALSE)
        box()}
    n.groups <- length(levels(groups))
    if (n.groups > length(col)) stop("number of groups exceeds number of available colors")
    if (length(col) == 1) col <- rep(col, 3)
    indices <- NULL
    range.x <- if (logged("x")) range(log(.x), na.rm=TRUE) else range(.x, na.rm=TRUE)
    counts <- table(groups)
    if (any(counts == 0)){
        levels <- levels(groups)
        warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
    }
    for (i in 1:n.groups){
        if (counts[i] == 0) next
        subs <- groups == levels(groups)[i]
        points(if (is.null(jitter$x) || jitter$x == 0) .x[subs] else jitter(.x[subs], factor=jitter$x), 
               if (is.null(jitter$y) || jitter$y == 0) .y[subs] else jitter(.y[subs], factor=jitter$y), 
               pch=pch[i], col=col[if (n.groups == 1) 3 else i], cex=cex)
        if (by.groups){
            if (is.function(smoother)) smoother(.x[subs], .y[subs], col=col[i], 
                                                log.x=logged("x"), log.y=logged("y"), spread=spread, smoother.args=smoother.args)
            if (is.function(reg.line)) reg(reg.line, .x[subs], .y[subs], lty=lty, lwd=lwd, log.x=logged("x"), log.y=logged("y"), col=col[i])
            if (ellipse) {
                X <- na.omit(data.frame(x=.x[subs], y=.y[subs]))
                if (logged("x")) X$x <- log(x)
                if (logged("y")) X$y <- log(y)
                with(X, dataEllipse(x, y, plot.points=FALSE, lwd=1, log=log,
                                    levels=levels, col=col[i], robust=robust))
            }
            if (id.method[1] != "identify")
                indices <- c(indices,
                     showLabels(.x[subs], .y[subs], labels=labels[subs], id.method=id.method,
                     id.n=id.n, id.cex=id.cex, id.col=col[i], all=list(labels=labels, subs=subs)))
#            if (id.method[1] != "identify") indices <- c(indices,
#                                                         showLabels(.x[subs], .y[subs], labels=labels[subs], id.method=id.method,
#                                                                    id.n=id.n, id.cex=id.cex, id.col=col[i]))
        }}
    if (!by.groups){
        if (is.function(smoother)) smoother(.x, .y, col=col[2], 
                                            log.x=logged("x"), log.y=logged("y"), spread, smoother.args=smoother.args)
        if (is.function(reg.line)) reg(reg.line, .x, .y, lty=lty, lwd=lwd, log.x=logged("x"), log.y=logged("y"), col=col[1])
        if (ellipse) {
            X <- na.omit(data.frame(x=.x, y=.y))
            if (logged("x")) X$x <- log(X$x)
            if (logged("y")) X$y <- log(X$y)
            with(X, dataEllipse(x, y, plot.points=FALSE, lwd=1, log=log, levels=levels, col=col[1],
                                robust=robust))
        }
        if (id.method[1] != "identify") indices <- showLabels(
            .x, .y, labels=labels, 
            id.method=id.method, id.n=id.n, id.cex=id.cex, id.col=id.col)
    }
    if (legend.plot) {
        xpd <- par(xpd=TRUE)
        on.exit(par(xpd=xpd), add=TRUE)
        usr <- par("usr")
        if (missing(legend.coords)){
            legend.x <- if (logged("x")) 10^(usr[1]) else usr[1]
            legend.y <- if (logged("y")) 10^(usr[4] + 1.2*top*strheight("x")) else usr[4] + 1.2*top*strheight("x")
            legend.coords <- list(x=legend.x, y=legend.y)
        }
        legend(legend.coords, legend=levels(groups)[counts > 0], 
               pch=pch[counts > 0], col=col[1:n.groups][counts > 0], pt.cex=cex, cex=cex.lab, 
               title=legend.title, bg="white", ncol=legend.columns)
    }
    if (id.method[1] == "identify") indices <- showLabels(.x, .y, labels, 
                                                          id.method=id.method, id.n=length(.x), id.cex=id.cex, id.col=id.col)
    if (is.null(indices)) invisible(indices) else if (is.numeric(indices)) sort(indices) else indices
} 

sp <- function(x, ...)  UseMethod("scatterplot", x)


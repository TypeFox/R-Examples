"plot.relimplm" <-
function (x, ..., names.abbrev = 4, ylim=NULL, main=NULL, cex.title=1.5) 
{
    # function shows barplots
    if (!(is(x, "relimplm"))) 
        stop("x must be the output from function calc.relimp")
    if (!(is.numeric(names.abbrev))) 
        stop("names.abbrev must be a number")
    p <- length(x@namen) - 1
    yname <- x@namen[1]
    if (nchar(yname)>60) yname <- "response variable"
    if (is.null(main)) main <- paste("Relative importances for ",  yname, sep="")
    xnames <- substr(x@namen[2:(p + 1)], 1, names.abbrev)
    if (x@rela) ylab <- expression("% of " * R^2)
    else ylab <- "% of response variance"
    if (x@rela && !is.null(x@always)) ylab <- expression("% of remaining " * R^2)
    if (is.null(ylim)){
    maxi <- max(x@lmg, x@pmvd, x@last, x@first, x@betasq, x@pratt, x@genizi, x@car)
    mini <- min(0, x@lmg, x@pmvd, x@last, x@first, x@betasq, x@pratt, x@genizi, x@car)
        if (maxi>=0.5) axmax <- 10 * ceiling(10 * maxi)
        else {
        if (maxi>=0.1) axmax <- 5 * ceiling(20 * maxi)
            else axmax <- 2 * ceiling(50 * maxi)}
        axmin <- 2 * floor(50 * mini)}
        else {axmax<-ylim[2]; axmin<-ylim[1]}
    type <- x@type
    reltext="%, metrics are not normalized."
    if (x@rela) reltext = "%, metrics are normalized to sum 100%."
    subtext = " "
    if (!is.null(x@always)) {
           if (length(x@always)>1) subtext <- 
               paste("Analysis adjusted for ", length(x@alwaysnam), " regressors that make up ", 
               round(100 * (x@R2 - x@R2.decomp),2), " pct. pts. of " , sep = "")
           if (length(x@always)==1) subtext <- 
               paste("Analysis adjusted for one regressor that makes up ", 
               round(100 * (x@R2 - x@R2.decomp),2), " pct. pts. of " , sep = "")
           }
    xlab <- " "

    if (length(type) == 0) 
        print("Nothing to plot")
    else {
        ntype <- length(type)
        op <- par(no.readonly = TRUE)
        oma <- c(3, 0, 3, 0)
        par(mar = c(2.6, 4.1, 4.1, 2.1))
        if (!is.null(x@always)) {
               oma[1] <- 4.5
               } 
        if (ntype == 1) 
            par(oma=oma)
        if (ntype == 2) 
            par(mfrow = c(1, 2), oma=oma)
        if (ntype > 2 && ntype <= 4) 
            par(mfrow = c(2, 2), oma=oma)
        if (ntype > 2 && ntype <= 4) 
            par(mfrow = c(2, 2), oma=oma)
        if (ntype > 4) 
            par(mfrow = c(2, 3), oma=oma)
        if (ntype > 6) 
            par(mfrow = c(2, 4), oma=oma)
        if ("lmg" %in% x@type) 
            barplot(100 * x@lmg, ylab = ylab, main = "Method LMG",
                names.arg = xnames, ylim = c(axmin, axmax), xlab = xlab)
        if ("pmvd" %in% x@type) 
            barplot(100 * x@pmvd, ylab = ylab, main = "Method PMVD",
                names.arg = xnames, ylim = c(axmin, axmax), xlab = xlab)
        if ("last" %in% x@type) 
            barplot(100 * x@last, ylab = ylab, main = "Method Last",
                names.arg = xnames, ylim = c(axmin, axmax), xlab = xlab)
        if ("first" %in% x@type) 
            barplot(100 * x@first, ylab = ylab, main = "Method First",
                names.arg = xnames, ylim = c(axmin, axmax), xlab = xlab)
        if ("betasq" %in% x@type) 
            barplot(100 * x@betasq, ylab = ylab, main = "Method Betasq", 
                names.arg = xnames, ylim = c(axmin, axmax), xlab = xlab)
        if ("pratt" %in% x@type) 
            barplot(100 * x@pratt, ylab = ylab, main = "Method Pratt", 
                names.arg = xnames, ylim = c(axmin, axmax), xlab = xlab)
        if ("genizi" %in% x@type) 
            barplot(100 * x@genizi, ylab = ylab, main = "Method Genizi", 
                names.arg = xnames, ylim = c(axmin, axmax), xlab = xlab)
        if ("car" %in% x@type) 
            barplot(100 * x@car, ylab = ylab, main = "Method CAR", 
                names.arg = xnames, ylim = c(axmin, axmax), xlab = xlab)
        title(main=main, outer=T, cex.main=cex.title)
        mtext(bquote(R^2==.(100*round(x@R2,4)) * .(eval(reltext))), side = 1, line=1.5, outer=T, cex=1, adj=0.5)
        if (!is.null(x@always)) 
            mtext(bquote(.(eval(subtext)) * R^2), side = 1, line=3, outer=T, cex=1, adj=0.5)
        par(op)
    }
}


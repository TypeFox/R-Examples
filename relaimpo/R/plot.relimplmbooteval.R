"plot.relimplmbooteval" <- 
function (x, ..., lev = max(x@level), names.abbrev = 4, ylim=NULL, main=NULL, cex.title=1.5) 
{
    # function shows barplots with error bars indicating confidence interval for chosen level
    # if chosen level not available, confidence interval for largest available level is produced

    #error control
    if (!(is(x, "relimplmbooteval") || is(x, "relimplmbootMI"))) 
        stop("x must be the output from function booteval.relimp or mianalyze.relimp")
    if (!(is.numeric(names.abbrev))) 
        stop("names.abbrev must be a number")
    if (!(is.numeric(lev))) 
        stop("lev must be a number")

    #no good way of tilting available
    #vertical labels do not work as desired but overlap with sub text
    #horizontal plotting does not work either
    #current solution: abbreviation of names to at most names.abbrev characters, default 4

    p <- length(x@namen) - 1
    yname <- x@namen[1]
    if (nchar(yname)>60) yname <- "response variable"
    if (is.null(main)) main <- paste("Relative importances for ",  yname, sep="")
    xnames <- substr(x@namen[2:(p + 1)], 1, names.abbrev)
    level <- x@level
    pick <- which(level == lev)
    if (length(pick) == 0) {
        pick <- which.max(level)
        cat("Chosen confidence level ", 100 * lev, "% not available,", 
            "\n")
        cat("largest available level (=default) plotted instead.", 
            "\n")
        cat("Available levels: ", level, "\n", sep = " ")
    }
    if (x@rela) ylab <- expression("% of " * R^2)
    else ylab <- "% of response variance"
    if (x@rela && !is.null(x@always)) ylab <- expression("% of remaining " * R^2)
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
    type <- x@type
    if (length(type) == 0) 
        print("Nothing to plot")
    else {
        if (is.null(ylim)){ 
        maxi <- 0
        mini <- 0
        for (a in type) {
            maxi <- max(maxi, slot(x, paste(a, "upper", sep = "."))[pick, 
                ])
            mini <- min(mini, slot(x, paste(a, "lower", sep = "."))[pick, 
                ])
        }
        if (maxi>=0.5) axmax <- 10 * ceiling(10 * maxi)
        else {
        if (maxi>=0.1) axmax <- 5 * ceiling(20 * maxi)
            else axmax <- 2 * ceiling(50 * maxi)}
        axmin <- 2 * floor(50 * mini)}
        else {axmax<-ylim[2]; axmin<-ylim[1]}
        ntype <- length(type)
        op <- par(no.readonly = TRUE)
        oma <- c(3, 0, 4, 0)
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

        if ("lmg" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@lmg, decreasing = T, index = T)$ix
            plt <- barplot(100 * x@lmg[index], main = "Method LMG",
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab, xlab = xlab)
            segments(plt, 100 * x@lmg.lower[pick, ][index], plt, 100 * x@lmg.upper[pick, 
                ][index])
        }
        if ("pmvd" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@pmvd, decreasing = T, index = T)$ix
            plt <- barplot(100 * x@pmvd[index], main = "Method PMVD", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab, xlab = xlab)
            segments(plt, 100 * x@pmvd.lower[pick, ][index], plt, 100 * x@pmvd.upper[pick, 
                ][index])
        }
        if ("last" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@last, decreasing = T, index = T)$ix
            plt <- barplot(100 * x@last[index], main = "Method Last", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab, xlab = xlab)
            segments(plt, 100 * x@last.lower[pick, ][index], plt, 100 * x@last.upper[pick, 
                ][index])
        }
        if ("first" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@first, decreasing = T, index = T)$ix
            plt <- barplot(100 * x@first[index], main = "Method First", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab, xlab = xlab)
            segments(plt, 100 * x@first.lower[pick, ][index], plt, 
                100 * x@first.upper[pick, ][index])
        }
        if ("betasq" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@betasq, decreasing = T, index = T)$ix
            plt <- barplot(100 * x@betasq[index], main = "Method Betasq", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab, xlab = xlab)
            segments(plt, 100 * x@betasq.lower[pick, ][index], plt, 
                100 * x@betasq.upper[pick, ][index])
        }
        if ("pratt" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@pratt, decreasing = T, index = T)$ix
            plt <- barplot(100 * x@pratt[index], main = "Method Pratt", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab, xlab = xlab)
            segments(plt, 100 * x@pratt.lower[pick, ][index], plt, 
                100 * x@pratt.upper[pick, ][index])
        }
        if ("genizi" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@genizi, decreasing = T, index = T)$ix
            plt <- barplot(100 * x@genizi[index], main = "Method Genizi", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab, xlab = xlab)
            segments(plt, 100 * x@genizi.lower[pick, ][index], plt, 
                100 * x@genizi.upper[pick, ][index])
        }
        if ("car" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@car, decreasing = T, index = T)$ix
            plt <- barplot(100 * x@car[index], main = "Method CAR", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab, xlab = xlab)
            segments(plt, 100 * x@car.lower[pick, ][index], plt, 
                100 * x@car.upper[pick, ][index])
        }
        title(main=main, line=2, outer=T, cex.main=cex.title)
        if (!x@fixed) title(main=paste("with ", 100 * level[pick], 
                  "% bootstrap confidence intervals", sep=""), line=0.5, outer=T, cex.main=1.2)
        if (x@fixed) title(main=paste("with ", 100 * level[pick], 
                  "% bootstrap confidence intervals (X-matrix treated as fixed)", sep=""), line=0.5, outer=T, 
                  cex.main=1.2)
        mtext(bquote(R^2==.(100*round(x@R2,4)) * .(eval(reltext))), side = 1, line=1.5, outer=T, cex=1, adj=0.5)
        if (!is.null(x@always)) 
            mtext(bquote(.(eval(subtext)) * R^2), side = 1, line=3, outer=T, cex=1, adj=0.5)
        par(op)
    }
}


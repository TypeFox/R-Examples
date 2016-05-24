# From Hmisc::wtd.var
wtd.sd <- function(x, weights)
{
    xbar <- sum(weights * x)/sum(weights)
    sqrt(sum(weights * ((x - xbar)^2))/(sum(weights) - 1))
}


compmeans <- function(x, f, w, sort = FALSE, maxlevels = 60,
                      user.missing, missing.include = FALSE,
                      plot = getOption("descr.plot"),
                      relative.widths = TRUE, col = "lightgray",
                      warn = getOption("descr.warn"), ...)
{
    row.label <- attr(f, "label")
    column.label <- attr(x, "label")
    row.name <- deparse(substitute(f))
    column.name <- deparse(substitute(x))

    f.name <- deparse(substitute(f))
    n.name <- deparse(substitute(x))
    lf <- length(f)
    lx <- length(x)
    if (lf != lx) {
        msg <- paste(f.name, gettext("and", domain = "R-descr"), n.name,
                     gettext("have different lengths", domain = "R-descr"))
        stop(msg)
    }

    if (is.factor(f) == FALSE) {
        f <- factor(f)
        nl <- length(levels(f))
        if (nl > maxlevels) {
            msg <- paste(f.name,
                         gettext("was converted into a factor, but the new variable had too many levels",
                                 domain = "R-descr"))
            stop(msg)
        }
        if(warn){
            wmsg <- paste(gettext("Warning:", domain = "R-descr"), " \"", f.name,
                          "\" ", gettext("was converted into factor!",
                                         domain = "R-descr"), sep = "")
            warning(wmsg)
        }
    } else{
        class(f) <- "factor"
    }
    if(!missing(user.missing)){
        user.missing <- paste("^", user.missing, "$", sep = "")
        flevels <- levels(f)
        for(lev in user.missing){
            if(length(grep(lev, flevels))){
                idx <- grep(lev, as.character(f)) 
                if(length(idx))
                    f[idx] <- NA
            }
        }
        f <- factor(f)
    }
    if(missing.include)
        f <- no.drop.levels(f)

    if(is.numeric(x)){
        class(x) <- "numeric"
    }
    if (missing(w)) {
        wt <- rep(1, lf)
    } else {
        wt <- w
        lw <- length(w)
        if (lw != lf) {
            msg <- paste(f.name, gettext("and", domain = "R-descr"), "weight",
                         gettext("have different lengths.", domain = "R-descr"))
            stop(msg)
        }
    }
    if(is.numeric(wt)){
        class(wt) <- "numeric"
    }

    if (is.factor(x) == TRUE) {
        x <- as.numeric(x)
        if(warn){
            wmsg <- paste(gettext("Warning:", domain = "R-descr"), " \"", n.name, "\" ", 
                          gettext("was converted from factor into numeric!", domain = "R-descr"))
            warning(wmsg)
        }
    }
    has.w <- FALSE
    k <- grep(FALSE, (is.na(f) | is.na(x) | is.na(wt)))
    f <- f[k]
    x <- x[k]
    wt <- wt[k]
    lf2 <- length(f)
    if (lf > lf2 && warn) {
        cat("\n")
        msg <- gettext("rows with missing values dropped", domain = "R-descr")
        wmsg <- paste(lf - lf2, msg)
        warning(wmsg)
    }

    xwsum <- tapply(x * wt, f, sum)
    wsum <- tapply(wt, f, sum)
    xmean <- xwsum / wsum
    wsum <- round(wsum)
    wsd <- xmean

    nflevs <- length(levels(f))
    b <- split(data.frame(x, wt), f)
    wsd <- sapply(b, function(.df) wtd.sd(.df$x, .df$wt))

    width <- wsum
    l <- length(xmean)
    xmean[l+1] <- weighted.mean(x, wt)
    wsum[l+1] <- round(sum(wt))
    wsd[l+1] <- wtd.sd(x, wt)
    tab <- cbind(xmean, wsum, wsd)
    tabrn <- rownames(tab)
    tabrn[l+1] <- gettext("Total", domain = "R-descr")
    rownames(tab) <- tabrn
    colnames(tab) <- c(gettext("Mean", domain = "R-descr"),
                       gettext("N", domain = "R-descr"),
                       gettext("Std. Dev.", domain = "R-descr"))
    if(sort){
        len <- length(xmean)
        len1 <- len - 1
        ordl <- order(xmean[1:len1]) # Do not sort the "Total"
        tab <- tab[c(ordl, len),]
        width <- width[ordl]
        f <- factor(as.numeric(f), levels = ordl, labels = levels(f)[ordl])
    }
    attr(tab, "row.name") <- row.name
    attr(tab, "column.name") <- column.name
    attr(tab, "row.label") <- row.label
    attr(tab, "column.label") <- column.label
    attr(tab, "maxlevels") <- maxlevels
    attr(tab, "col") <- col
    if(relative.widths)
        attr(tab, "width") <- width
    else
        attr(tab, "width") <- rep(1, length(width))

    # Add attributes to plot the object:
    attr(tab, "x") <- x
    attr(tab, "f") <- f
    if(!missing(w))
        attr(tab, "wt") <- wt
    class(tab) <- c("meanscomp", "matrix")

    if(plot)
        plot.meanscomp(tab, ...)
    tab
}

print.meanscomp <- function(x, ...)
{
    rlab <- attr(x, "row.label")
    clab <- attr(x, "column.label")

    if(is.null(rlab))
        rlab <- attr(x, "row.name")
    if(is.null(clab))
        clab <- attr(x, "column.name")

    # 'domain' is necessary because this function is not exported to
    # 'descr' namespace.
    msg1 <- gettext("Mean value of", domain = "R-descr")
    msg2 <- gettext("according to", domain = "R-descr")
    lwd <- getOption("width")
    msg <- paste(msg1, ' "', clab, '" ', msg2, ' "', rlab, '"', sep = "")

    # Break the label string if it is too large:
    if(nchar(msg, type = "width") < lwd){
        cat(msg, "\n", sep = "")
    } else {
        if((nchar(msg1, type = "width") + nchar(clab, type = "width")) < lwd) {
            msg <- paste(msg1, ' "', clab, '" ', sep = "")
            if((nchar(msg, type = "width") + nchar(msg2, type = "width")) < lwd) {
                cat(msg, msg2, '\n', '"', rlab, '"', '\n', sep = "")
            } else {
                cat(msg, "\n", sep = "")
                if((nchar(msg2, type = "width") + nchar(rlab, type = "width")) < (lwd - 1)){
                    cat(msg2, ' "', rlab, '"\n', sep = "")
                } else {
                    cat(msg2, '\n"', rlab, '"\n', sep = "")
                }
            }
        } else {
            cat(msg1, '\n"', clab, '"\n', sep = "")
            if((nchar(msg2, type = "width") + nchar(rlab, type = "width")) < (lwd - 1)){
                cat(msg2, ' "', rlab, '"\n', sep = "")
            } else {
                cat(msg2, '\n"', rlab, '"\n', sep = "")
            }
        }
    }
    attr(x, "row.name") <- NULL
    attr(x, "column.name") <- NULL
    attr(x, "row.label") <- NULL
    attr(x, "column.label") <- NULL
    attr(x, "maxlevels") <- NULL
    attr(x, "x") <- NULL
    attr(x, "f") <- NULL
    attr(x, "wt") <- NULL
    attr(x, "col") <- NULL
    attr(x, "width") <- NULL
    class(x) <- "matrix"
    print(x, ...)
    return(invisible(NULL))
}

plot.meanscomp <- function(x, xlab, ylab, width, col, ...)
{
    if(missing(xlab))
        xlab <- attr(x, "row.name")
    if(missing(ylab))
        ylab <- attr(x, "column.name")
    if(missing(width))
        width <- attr(x, "width")
    if(missing(col))
        col <- attr(x, "col")
    maxlevels <- attr(x, "maxlevels")
    v <- attr(x, "x")
    f <- attr(x, "f")
    w <- attr(x, "wt")

    if(!is.factor(f))
        stop(gettext("f is not a factor.", domain = "R-descr"))
    if(length(levels(f)) > maxlevels)
        stop(gettext("Number of levels of \"f\" is higher than maxlevels.",
                     domain = "R-descr"))

    if(is.null(w)){
        boxplot(v ~ f, ylab = ylab, xlab = xlab, width = width, col = col, ...)
    } else {
        d <- data.frame(v, f, w)
        dd <- split(d, f)
        zz <- lapply(dd, function(d) spBwplotStats(d$v, d$w))
        z <- zz[[1]]
        z$nzero <- NULL
        z <- unclass(z)
        z$group <- rep(1, length(z$out))
        z$names <- levels(f)
        for(i in 2:length(zz)){
            z$stats <- cbind(z$stats, zz[[i]]$stats)
            z$n <- c(z$n, zz[[i]]$n)
            z$out <- c(z$out, zz[[i]]$out)
            z$group <- c(z$group, rep(i, length(zz[[i]]$out)))
        }
        bxp(z, ylab = ylab, xlab = xlab, width = width, boxfill = col, ...)
        return(invisible(z))
    }
}

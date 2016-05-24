
freq <- function (x, w, user.missing, plot = getOption("descr.plot"), ...)
{
    xlab <- attr(x, "label", TRUE)
    if(is.null(xlab))
        xlab <- deparse(substitute(x))
    if (is.factor(x) == FALSE)
        x <- as.factor(x)
    xclass <- class(x)

    if (missing(w))
        w <- rep(1, length(x))

    nmiss <- sum(is.na(x))
    xlevels <- levels(x)
    l <- length(xlevels)
    hasna <- FALSE
    xv <- x
    if (nmiss) {
        hasna <- TRUE
        l <- l + 1
        xlevels[l] <- "NA's"
        x <- as.numeric(x)
        x[is.na(x)] <- l
        x <- factor(x, levels=1:l, labels = xlevels)
    }

    xfreq <- tapply(w, x, sum, na.rm = TRUE)
    xfreq[is.na(xfreq)] <- 0
    xtotal <- sum(xfreq, na.rm = TRUE)
    xperc <- 100 * xfreq / xtotal

    ftab <- cbind(xfreq, xperc)
    cnames <- c(gettext("Frequency", domain = "R-descr"),
                gettext("Percent", domain = "R-descr"))

    xvfreq <- xfreq
    if(nmiss){
        xvfreq[xlevels == "NA's"] <- NA
    }
    if(!missing(user.missing)){
        user.missing <- paste("^", user.missing, "$", sep = "")
        for(lev in user.missing){
            idx <- grep(lev, xlevels)
            if(length(idx))
                xvfreq[idx] <- NA
        }
    }

    if(nmiss || !missing(user.missing)){
        xvtotal <- sum(xvfreq, na.rm = TRUE)
        xvperc <- 100 * xvfreq / xvtotal
        ftab <- cbind(ftab, xvperc)
        cnames <- c(cnames, gettext("Valid Percent", domain = "R-descr"))
    }

    if(xclass[1] == "ordered"){
        if(nmiss || !missing(user.missing)){
            xxvperc <- xvperc
            xxvperc[is.na(xxvperc)] <- 0
            xvcumsum <- cumsum(xxvperc)
            xvcumsum[is.na(xvperc)] <- NA
        } else
            xvcumsum <- cumsum(xperc)
        ftab <- cbind(ftab, xvcumsum)
        cnames <- c(cnames, gettext("Cum Percent", domain = "R-descr"))
    }

    total <- apply(ftab, 2, sum, na.rm = TRUE)
    if(xclass[1] == "ordered")
        total["xvcumsum"] <- NA
    ftab <- rbind(ftab, total)

    rnames <- levels(x)
    rnames[l + 1] <- gettext("Total", domain = "R-descr")

    colnames(ftab) <- cnames
    rownames(ftab) <- rnames

    attr(ftab, "xlab") <- xlab
    class(ftab) <- c("freqtable", "matrix")

    # Attributes for plotting
    if(nmiss || !missing(user.missing))
        xdata.c <- xvfreq
    else
        xdata.c <- xfreq
    if(length(grep("^NA's$", names(xdata.c))) > 0)
        xdata.c["NA's"] <- NA
    xdata.c <- xdata.c[!is.na(xdata.c)]
    if(nmiss || !missing(user.missing))
        xdata.p <- xvperc
    else
        xdata.p <- xperc
    if(length(grep("^NA's$", names(xdata.p))) > 0)
        xdata.p["NA's"] <- NA
    xdata.p <- xdata.p[!is.na(xdata.p)]

    attr(ftab, "xdata.c") <- xdata.c
    attr(ftab, "xdata.p") <- xdata.p

    if(plot == TRUE)
        plot.freqtable(ftab, ...)

    ftab
}

print.freqtable <- function(x, digits = 4, na.print="", ...){
    xlab <- attr(x, "xlab")
    cat(xlab, "\n")
    attr(x, "xlab") <- NULL
    attr(x, "xdata.c") <- NULL
    attr(x, "xdata.p") <- NULL
    class(x) <- "matrix"
    print(x, digits = digits, na.print = na.print, ...)
    return(invisible(NULL))
}

plot.freqtable <- function(x, y.axis = "count", ...)
{
    if(y.axis == "count"){
        xdata <- attr(x, "xdata.c")
    } else if(y.axis == "percent"){
        xdata <- attr(x, "xdata.p")
    } else {
        msg <- paste(gettext("Invalid y.axis: '", domain = "R-descr"),
                     y.axis[1], "'", sep = "")
        stop(msg)
    }
    barplot(xdata, ...)
}


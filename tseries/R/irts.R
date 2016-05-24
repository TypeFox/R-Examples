## Copyright (C) 1997-2003  Adrian Trapletti
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

##
## Irregular time series objects
##

value <- function (x, ...) UseMethod ("value")
                    
irts <-
function(time, value)
{
    if(inherits(time, "POSIXct")) {
        time <- as.numeric(time)
    }
    if(!is.vector(time))
        stop("time is not a vector")
    if(!is.vector(value) && !is.matrix(value))
        stop("value is not a vector and not a matrix")
    if(length(time) != NROW(value))
        stop("time and value have not the same number of rows")
    class(time) <- c("POSIXt", "POSIXct")
    irts <- list(time = time, value = value)
    class(irts) <- "irts"
    return(irts)
}

is.irts <-
function(object)
{
    return(inherits(object, "irts"))
}

as.irts <- function(object) UseMethod("as.irts")

as.irts.default <-
function(object)
{
    return(irts(object[,1], object[,-1]))
}

as.irts.zoo <-
function(object, ...)
{
    index <- attr(object, "index")
    stopifnot(inherits(index, "POSIXct"))
    attr(object, "index") <- NULL
    irts(index, unclass(object))  
}

value.irts <-
function(x, ...)
{
    if(!inherits(x, "irts"))
        stop("method is only for irts objects")
    return(x$value)
}

time.irts <-
function(x, ...)
{
    if(!inherits(x, "irts"))
        stop("method is only for irts objects")
    return(x$time)
}

print.irts <-
function(x, format = "%Y-%m-%d %H:%M:%S", tz = "GMT",
         usetz = TRUE, format.value = NULL, ...)
{
    if(!inherits(x, "irts"))
        stop("method is only for irts objects")
    n <- length(x$time)
    for(i in 1:n) {
        cat(format(x$time[i], format = format, tz = tz, usetz = usetz))
        cat(" ")
        if(is.vector(x$value))
            cat(formatC(x$value[i], format = format.value, ...))
        else
            cat(formatC(x$value[i,], format = format.value, ...))
        cat("\n")
    }
    invisible(x)
}

read.irts <-
function(file, format = "%Y-%m-%d %H:%M:%S", tz = "GMT", ...)
{
    seqN <- function(from, to) {
        if((0 == length(from)) || (0 == length(to)))
            NULL
        else if(to-from+1 <= 0) 
            NULL
        else seq(from, to)
    }
    
    data <- read.table(file, as.is = TRUE, ...)
    n <- length(unlist(strsplit(format, split = " ")))
    tmp <- data[,1]
    j <- 2
    while(j <= n) {
        tmp <- paste(tmp, data[,j])
        j <- j+1
    }
    time <- as.numeric(as.POSIXct(strptime(tmp, format = format), tz = tz))
    value <- as.matrix(data[,-seqN(1, n)])
    return(irts(time, value[,,drop = TRUE]))
}

write.irts <-
function(object, file = "", append = FALSE, quote = FALSE, sep = " ", eol = "\n",
         na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = "escape",
         format = "%Y-%m-%d %H:%M:%S", tz = "GMT", usetz = FALSE, format.value = NULL, ...)
{
    dataframe <- data.frame(time = format(object$time, format = format, tz = tz, usetz = usetz),
                            value = formatC(object$value, format = format.value, ...))
    write.table(dataframe, file = file, append = append, quote = quote, sep = sep, eol = eol,
                na = na, dec = dec, row.names = row.names, col.names = col.names, qmethod = qmethod)
    invisible(object)
}

weekday <-
function(object, tz = "GMT")
{
    if(!inherits(object, "irts"))
        stop("function is only for irts objects")
    return(as.POSIXlt(object$time, tz = tz)$wday)
}

daysecond <-
function(object, tz = "GMT")
{
    if(!inherits(object, "irts"))
        stop("function is only for irts objects")
    hour <- as.POSIXlt(object$time, tz = tz)$hour
    min <- as.POSIXlt(object$time, tz = tz)$min
    sec <- as.POSIXlt(object$time, tz = tz)$sec
    return(3600*hour+60*min+sec)
}

is.businessday <-
function(object, tz = "GMT")
{
    if(!inherits(object, "irts"))
        stop("function is only for irts objects")
    wday <- as.POSIXlt(object$time, tz = tz)$wday
    return((0 < wday) & (wday < 6))
}

is.weekend <- function(object, tz = "GMT")
{
    if(!inherits(object, "irts"))
        stop("function is only for irts objects")
    wday <- as.POSIXlt(object$time, tz = tz)$wday
    return((0 == wday) | (wday == 6))
}

"[.irts" <-
function(x, i, j, ...)
{
    if(!inherits(x, "irts"))
        stop("method is only for irts objects")
    if(is.vector(x$value)) {
        if(nargs() > 2) {   
            stop("incorrect number of dimensions")
        }
        if(missing(i)) {
            return(x)
        } else {
            return(irts(as.numeric(x$time)[i], x$value[i]))
        }
    } else {
        if(missing(i)) {
            if(missing(j)) {
                return(x)
            } else {
                return(irts(as.numeric(x$time), x$value[,j,drop = FALSE]))
            }
        } else {
            if(missing(j)) {
                return(irts(as.numeric(x$time)[i], x$value[i,,drop = FALSE]))
            } else {
                return(irts(as.numeric(x$time)[i], x$value[i,j,drop = FALSE]))
            }
        }
    }
}

approx.irts <-
function(object, time, ...)
{
    if(!inherits(object, "irts"))
        stop("function is only for irts objects")
    if(!inherits(time, "POSIXct"))
        stop("time is not of class POSIXct")
    ovalue <- as.matrix(object$value)
    otime <- as.numeric(object$time)
    time <- as.numeric(time)
    value <- matrix(0, NROW(time), NCOL(ovalue))
    for(i in 1:NCOL(ovalue)) {
        result <- approx(otime, ovalue[,i,drop = TRUE], time, ...)
        value[,i] <- result$y
    }
    return(irts(time, value[,,drop = TRUE]))
}

plot.irts <-
function(x, type = "l", plot.type = c("multiple", "single"),
         xlab = "Time", ylab = NULL, main = NULL, ylim = NULL,
         oma = c(6, 0, 5, 0), ...)
{
    seqN <- function(from, to) {
        if((0 == length(from)) || (0 == length(to)))
            NULL
        else if(to-from+1 <= 0) 
            NULL
        else seq(from, to)
    }
    
    addmain <- function(main, cex.main = par("cex.main"),
                        font.main = par("font.main"), 
                        col.main = par("col.main"), ...) {
        mtext(main, 3, 3, cex = cex.main, font = font.main, col = col.main, ...)
    }
    
    if(!inherits(x, "irts"))
        stop("method is only for irts objects")
    t <- time(x)
    v <- value(x)
    nser <- NCOL(v)
    if(is.null(main)) 
        main <- deparse(substitute(x))
    if(nser == 1) {
        if(is.null(ylab))
            ylab <- "Series"
        if(is.null(ylim))
            ylim <- range(v[is.finite(v)])
        plot(t, v, type = type, xlab = xlab, ylab = ylab,
             main = main, ylim = ylim, ...)
    } else if(nser <= 10) {
        plot.type <- match.arg(plot.type)
        if(is.null(ylab)) {
            ylab <- colnames(v)
            if(is.null(ylab)) 
                ylab <- paste("Series", 1:nser)
        }
        if(plot.type == "single") {
            if(is.null(ylim))
                ylim <- range(v[is.finite(v)])
            plot.default(t, v[,1], type = type, xlab = xlab, ylab = ylab,
                         main = main, ylim = ylim, xaxt = "n", ...) 
            for(i in seqN(2, nser)) {
                points(t, v[,i], type = type, xaxt = "n") 
            }
            axis.POSIXct(1, t)
        } else if(plot.type == "multiple") {
            oldpar <- par("mar", "oma", "mfcol")
            on.exit(par(oldpar))
            par(mar = c(0, 5.1, 0, 2.1), oma = oma)
            nc <- if(nser > 4) 2 else 1
            nr <- ceiling(nser/nc)
            par(mfcol = c(nr, nc))
            for(i in seqN(1, nser)) {
                plot.default(t, v[,i], type = type, xlab = xlab, ylab = "", xaxt = "n", ...)
                mtext(ylab[i], 2, 3)
                if((i%%nr == 0) || (i == nser))
                    axis.POSIXct(1, t)
            }
            if(!is.null(main)) {
                par(mfcol = c(1, 1))
                addmain(main, ...)
            }
        }
    } else {
        stop("cannot plot more than 10 series")
    }
    invisible(x)
}

lines.irts <-
function(x, type = "l", ...)
{
    if(!inherits(x, "irts"))
        stop("method is only for irts objects")
    t <- time(x)
    v <- value(x)
    nser <- NCOL(v)
    if(nser == 1) {
        lines(t, v, type = type, ...)
    } else {
        stop("cannot plot multivariate irregular time-series object")
    }
    invisible(x)
}

points.irts <-
function(x, type = "p", ...)
{
    if(!inherits(x, "irts"))
        stop("method is only for irts objects")
    t <- time(x)
    v <- value(x)
    nser <- NCOL(v)
    if(nser == 1) {
        points(t, v, type = type, ...)
    } else {
        stop("cannot plot multivariate irregular time-series object")
    }
    invisible(x)
}


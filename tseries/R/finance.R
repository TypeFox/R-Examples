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
## Financial time series analysis
##

portfolio.optim <- function (x, ...) UseMethod ("portfolio.optim")

portfolio.optim.ts <-
function (x, ...)
{
    if(!is.ts(x))
        stop("method is only for time series")
    if(NCOL(x) == 1)
        stop("x is not a multivariate time series")
    res <- portfolio.optim.default(as.matrix(x), ...)
    res$px <- ts(res$px, start = start(x), frequency = frequency(x))
    return(res)
}

portfolio.optim.default <-
function(x, pm = mean(x), riskless = FALSE, shorts = FALSE,
         rf = 0.0, reslow = NULL, reshigh = NULL, covmat = cov(x), ...) 
{
    if(NCOL(x) == 1)
        stop("x is not a matrix")
    if(any(is.na(x)))
        stop("NAs in x")
    k <- dim(x)[2]
    if(!is.matrix(covmat)) {
        stop("covmat is not a matrix")
    }
    if((dim(covmat)[1] !=k) || (dim(covmat)[2] !=k)) {
      stop("covmat has not the right dimension")
    }
    Dmat <- covmat
    dvec <- rep.int(0, k)
    big <- 1e+100
    if(!is.null(reslow) && is.null(reshigh)) {
        reshigh <- rep.int(big, k)
    }
    if(is.null(reslow) && !is.null(reshigh)) {
        reslow <- -rep.int(big, k)
    }
    if(!is.null(reslow)) {
        if(!is.vector(reslow)) {
            stop("reslow is not a vector")
        }
        if(length(reslow) != k) {
            stop("reslow has not the right dimension")
        }
    }
    if(!is.null(reshigh)) {
        if(!is.vector(reshigh)) {
            stop("reshigh is not a vector")
        }
        if(length(reshigh) != k) {
            stop("reshigh has not the right dimension")
        }
    }
    if(riskless) {
        a1 <- colMeans(x) - rf
        if(shorts) {
            a2 <- NULL
            b2 <- NULL
        }
        else {
            a2 <- matrix(0, k, k)
            diag(a2) <- 1
            b2 <- rep.int(0, k)
        }
        if(!is.null(reslow) && !is.null(reshigh)) {
            a3 <- matrix(0, k, k)
            diag(a3) <- 1
            Amat <- t(rbind(a1, a2, a3, -a3))
            b0 <- c(pm-rf, b2, reslow, -reshigh)
        }
        else {
            Amat <- t(rbind(a1, a2))
            b0 <- c(pm-rf, b2)
        }
        res <- solve.QP(Dmat, dvec, Amat, bvec=b0, meq=1)
    }
    else {
        a1 <- rep.int(1, k)
        a2 <- colMeans(x)
        if(shorts) {
            if(!is.null(reslow) && !is.null(reshigh)) {
                a3 <- matrix(0, k, k)
                diag(a3) <- 1
                Amat <- t(rbind(a1, a2, a3, -a3))
                b0 <- c(1, pm, reslow, -reshigh)
            }
            else {
                Amat <- t(rbind(a1, a2))
                b0 <- c(1, pm)
            }              
        }
        else {
            a3 <- matrix(0, k, k)
            diag(a3) <- 1
            b3 <- rep.int(0, k)
            if(!is.null(reslow) && !is.null(reshigh)) {
                Amat <- t(rbind(a1, a2, a3, a3, -a3))
                b0 <- c(1, pm, b3, reslow, -reshigh)
            }
            else {
                Amat <- t(rbind(a1, a2, a3))
                b0 <- c(1, pm, b3)
            }
        }
        res <- solve.QP(Dmat, dvec, Amat, bvec=b0, meq=2)
    }
    y <- c(tcrossprod(res$solution, x))
    ans <- list(pw = res$solution, px = y, pm = mean(y), ps = sd(y))
    return(ans)
}

get.hist.quote <-
function (instrument = "^gdax", start, end,
          quote = c("Open", "High", "Low", "Close"),
          provider = c("yahoo", "oanda"), method = NULL,
          origin = "1899-12-30", compression = "d",
	  retclass = c("zoo", "its", "ts"),
	  quiet = FALSE, drop = FALSE)
    ## Added new argument 'compression'.
    ## May be "d", "w" or "m", for daily weekly or monthly.
    ## Defaults to "d".
    ## John Bollinger, 2004-10-27, www.BollingerBands.com, bbands@yahoo.com
    ##
    ## Changed POSIXct class to Date class, 2005-03-31
    ## Added drop= arg, 2006-06-04
{
    if(missing(start)) start <- "1991-01-02"
    if(missing(end)) end <- format(Sys.Date() - 1, "%Y-%m-%d")
  
    provider <- match.arg(provider)
    retclass <- match.arg(retclass)

    start <- as.Date(start)
    end <- as.Date(end)

    if(is.null(method)) {
        method <- getOption("download.file.method")
        if(is.null(method))
            method <- "auto"
    }

    if(provider == "yahoo") {
        url <-
            paste("http://chart.yahoo.com/table.csv?s=",
                  instrument, 
                  format(start,
                         paste("&a=",
                               as.character(as.numeric(format(start, "%m"))-1),
                               "&b=%d&c=%Y",
                               sep = "")),
                  format(end,
                         paste("&d=",
                               as.character(as.numeric(format(end, "%m"))-1),
                               "&e=%d&f=%Y",
                               sep = "")), 
                  "&g=", compression,
                  "&q=q&y=0&z=", instrument,
                  "&x=.csv",
                  sep = "")
        destfile <- tempfile()

        ## Try downloading several times.
        ## Could make this controllable ...
        i <- 1L
        repeat {
            status <- tryCatch(download.file(url, destfile,
                                             method = method,
                                             quiet = quiet),
                               error = identity)
            if(!inherits(status, "error") && (status == 0))
                break
            unlink(destfile)            
            if(i >= 5L) {
                if(inherits(status, "error"))
                    stop(conditionMessage(status))
                else
                    stop(sprintf("download error, status %d", status))
            }
            message("download error, retrying ...")
            i <- i + 1L
        }
        
        nlines <- length(count.fields(destfile, sep = "\n"))
        if(nlines == 1) {
            unlink(destfile)
            stop(paste("no data available for", instrument))
        }
        
        ## Yahoo includes rows concerning dividends,
        ## hence need fill = TRUE and na.omit
        x <- read.table(destfile, header = TRUE, sep = ",", as.is = TRUE, fill = TRUE)
        x <- na.omit(x)

        ## As of 2012-07-07,
        ##   get.hist.quote("ibm",  start = "2012-07-06", end = "2012-07-06")
        ## gives two entries for 2012-07-06 ...
        ## Fix (hopefully only needed temporarily?) suggested by AZ.
        if(nrow(x) >= 2L && x[1L, 1L] == x[2L, 1L]) {
            warning("first date duplicated, first instance omitted")
            x <- x[-1L, , drop = FALSE]
        }
        
        ## Debug
        ## cat("read.table: start =", x[NROW(x),"Date"], "\n")
        ## cat("read.table: end   =", x[1,"Date"], "\n")
        
        unlink(destfile)
        
        names(x) <- gsub("\\.", "", names(x))
        nser <- pmatch(quote, names(x)[-1]) + 1
        if(any(is.na(nser)))
            stop("this quote is not available")
        n <- nrow(x)

        ## Yahoo now seems to format dates as %Y-%m-%d.
        dat <- as.Date(as.character(x[, 1]), "%Y-%m-%d")
        if(!quiet && (dat[n] != start))
            cat(format(dat[n], "time series starts %Y-%m-%d\n"))
        if(!quiet && (dat[1] != end))
            cat(format(dat[1], "time series ends   %Y-%m-%d\n"))

	if(retclass == "ts") {
            jdat <- unclass(julian(dat, origin = as.Date(origin)))
            ## We need unclass() because 1.7.0 does not allow adding a
            ## number to a "difftime" object.
            ind <- jdat - jdat[n] + 1
            y <- matrix(NA, nrow = max(ind), ncol = length(nser))
            y[ind, ] <- as.matrix(x[, nser, drop = FALSE])
            colnames(y) <- names(x)[nser]
	    y <- y[, seq_along(nser), drop = drop]
            return(ts(y, start = jdat[n], end = jdat[1]))
	} else {
	  x <- as.matrix(x[, nser, drop = FALSE])
	  rownames(x) <- NULL
	  y <- zoo(x, dat)
	  y <- y[, seq_along(nser), drop = drop]
	  if(retclass == "its") {
              if(inherits(tryCatch(getNamespace("its"), error = identity),
                          "error"))
                  warning("package its could not be loaded: zoo series returned")
              else {
                  index(y) <- as.POSIXct(index(y))
                  y <- its::as.its(y)
              }
	  }
	  return(y)
	}
    }
    else if(provider == "oanda") {
        if(!missing(quote)) {
            warning("argument 'quote' ignored for provider 'oanda'")
        }
        if(!missing(compression)) {
            warning("argument 'compression' ignored for provider 'oanda'")
        }

        currencies <- unlist(strsplit(instrument, split = "/"))
        ranges <- c("d7", "d30", "d60", "d90", "d180", "y1", "y2", "y5")
        range <- ranges[c(7, 30, 60, 90, 180, 364, 728, 1820) >=
                            difftime(end, start, units="days")][1]
        url <- paste0("https://www.oanda.com/currency/historical-rates/download?",
                      "quote_currency=", currencies[1],
                      "&end_date=", end,
                      "&start_date=", start,
                      "&period=daily",
                      "&display=absolute",
                      "&rate=0",
                      "&data_range=", range,
                      "&price=mid",
                      "&view=table",
                      "&base_currency_0=", currencies[2],
                      "&base_currency_1=",
                      "&base_currency_2=",
                      "&base_currency_3=",
                      "&base_currency_4=",
                      "&download=csv")
        destfile <- tempfile()
        
        status <- download.file(url, destfile, method = method, quiet = quiet)
        if(status != 0) {
            unlink(destfile)
            stop(paste("download error, status", status))
        }

        x <- read.csv(destfile, skip = 4, as.is = TRUE, header = TRUE)
        unlink(destfile)
        x <- head(x, -3)
        
        dat <- rev(as.Date(x[[1]]))
        n <- length(dat)
        if(!quiet && (dat[1] != start))
            cat(format(dat[1], "time series starts %Y-%m-%d\n"))
        if(!quiet && (dat[n] != end))
            cat(format(dat[n], "time series ends   %Y-%m-%d\n"))
        val <- rev(x[[2]])
        if(is.character(val))
            val <- as.numeric(sub(",", "", val, fixed = TRUE))

	if(retclass == "ts") {
            jdat <- unclass(julian(dat, origin = as.Date(origin)))
            ind <- jdat - jdat[1] + 1
            y <- rep.int(NA, max(ind))
            y[ind] <- val
            return(ts(y, start = jdat[1], end = jdat[n]))
	} else {
	  y <- zoo(val, dat)
	  if(retclass == "its") {
              if(inherits(tryCatch(getNamespace("its"), error = identity),
                          "error"))
                  warning("package its could not be loaded: zoo series returned")
              else {
                  index(y) <- as.POSIXct(index(y))
                  y <- its::as.its(y)
              }
	  }
	  return(y)
	}
    }
    else stop("provider not implemented")
}

maxdrawdown <-
function(x)
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    cmaxx <- cummax(x)-x
    mdd <- max(cmaxx)
    to <- which(mdd == cmaxx)
    from <- double(NROW(to))
    for (i in 1:NROW(to))
        from[i] <- max(which(cmaxx[1:to[i]] == 0))
    return(list(maxdrawdown = mdd, from = from, to = to))
}

sharpe <-
function(x, r = 0, scale = sqrt(250))
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(NROW(x) == 1)
        return(NA)
    else {
        y <- diff(x)
        return(scale * (mean(y)-r)/sd(y))
    }
}

sterling <-
function(x)
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(NROW(x) == 1)
        return(NA)
    else {
        return((x[NROW(x)]-x[1]) / maxdrawdown(x)$maxdrawdown)
    }
}

plotOHLC <-
function(x, xlim = NULL, ylim = NULL, xlab = "Time", ylab,
         col = par("col"), bg = par("bg"), axes = TRUE,
         frame.plot = axes, ann = par("ann"), main = NULL,
         date = c("calendar", "julian"), format = "%Y-%m-%d",
         origin = "1899-12-30", ...)
{
  if ((!is.mts(x)) ||
      (colnames(x)[1] != "Open") ||
      (colnames(x)[2] != "High") ||
      (colnames(x)[3] != "Low") ||
      (colnames(x)[4] != "Close"))
      stop("x is not a open/high/low/close time series")
  xlabel <- if (!missing(x)) 
      deparse(substitute(x))
  else NULL
  if (missing(ylab)) 
      ylab <- xlabel
  date <- match.arg(date)
  time.x <- time(x)
  dt <- min(lag(time.x)-time.x)/3
  if (is.null(xlim)) 
      xlim <- range(time.x)
  if (is.null(ylim)) 
      ylim <- range(x[is.finite(x)])
  plot.new()
  plot.window(xlim, ylim, ...)
  segments(time.x, x[, "High"], time.x, x[, "Low"], col = col[1], bg = bg)
  segments(time.x - dt, x[, "Open"], time.x, x[, "Open"],
           col = col[1], bg = bg)
  segments(time.x, x[, "Close"], time.x + dt, x[, "Close"],
           col = col[1], bg = bg)
  if (ann) 
      title(main = main, xlab = xlab, ylab = ylab, ...)  
  if (axes) {
      if (date == "julian") {
          axis(1, ...)
          axis(2, ...)
      }
      else {
          n <- NROW(x)
          lab.ind <- round(seq(1, n, length=5))
          D <- as.vector(time.x[lab.ind]*86400) + as.POSIXct(origin, tz = "GMT")
          DD <- format.POSIXct(D, format = format, tz ="GMT")
          axis(1, at=time.x[lab.ind], labels=DD, ...)
          axis(2, ...)
      }
  }
  if (frame.plot) 
      box(...)
}


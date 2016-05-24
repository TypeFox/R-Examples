mdcc <- function(chrono, clim, method = "response", start = 4, end =
                 9, timespan = NULL, vnames = NULL, sb = TRUE,
                 win.size = 25, win.offset = 1, startlast = TRUE,
                 boot = FALSE, ci = 0.05) {

  month.ids <- c(-1:-12, 1:12)
  errormsg1 <-
    "start and end have to define an interval in [-1, -2, ..., -12, 1, 2, ..., 12]."
  if (!is.element(start, month.ids) || !is.element(end, month.ids) ||
      which(month.ids == start) > which(month.ids == end)) {
                                        # check start and end dates
    stop(errormsg1)
  }
  ## climate data dispatcher gets called...
  clim <- climdispatch(clim)                    # properly formatted
                                        # climate data gets returned
                                        # here

  if (start*end > 0) {
    no.params <- (dim(clim)[2] -2)*length(start:end)
  } else {
    no.params <- (dim(clim)[2] -2)*length(start:end)-1 # 0 is not counted
  }
  ## raise error, when window size is smaller than number of params
  if (no.params > win.size) {
    win.size.msg <-
      paste("Window size is smaller than number of parameters! Consider adapting win.size to a minimum of ",
            no.params, ".", sep = "")
    stop(win.size.msg)
  }
  chrono.years <- as.numeric(row.names(chrono)) # get timespan of
                                        # chrono
  clim.years <- sort(unique(clim[, 1])) # get timespan of climate data
  if (chrono.years[1] <= clim.years[1]) {
    overlap <- na.omit(clim.years[match(chrono.years, clim.years)]) # get overlap of timespans
  } else {
    overlap <- na.omit(chrono.years[match(clim.years, chrono.years)])
  }
  if (is.null(timespan)) {
    start.year <- overlap[1]
    end.year <- tail(overlap, 1)
  } else {
    if (start > 0) {
      if (!is.element(timespan[1], overlap) ||
          !is.element(timespan[2], overlap)) {
        errormsg3 <- paste("timespan has to be between ", overlap[1],
                           " and ", tail(overlap, 1),
                           " for start dates in current year.",
                           sep = "")
        stop(errormsg3)
      } else {
        start.year <- timespan[1]
        end.year <- timespan[2]
      }
    } else {
      if (!is.element(timespan[1], overlap) ||
          !is.element(timespan[2], overlap)) {
        errormsg4 <- paste("timespan has to be between ", overlap[1] +
                           1, " and ", tail(overlap, 1),
                           " for start dates in previous year.", sep = "")
        stop(errormsg4)
      } else {
        start.year <- timespan[1]
        end.year <- timespan[2]
      }
    }
  }
  if (length(start.year:end.year) <= win.size) {
    errormsg5 <-
      paste("timespan is shorter than win.size. Consider adapting timespan to at least ",
            win.size, ".", sep = "")
    stop(errormsg5)
  }
  if (start < 0 && is.na(match((start.year - 1), clim.years))) { # check if a previous year is available in climatic data; otherwise set start.year + 1
    offset <- 1
  } else {
    offset <- 0
  }
  if (start < 0) { # make sure that data get truncated properly	
    interval.clim <-(start.year - 1 + offset):end.year
    interval.chrono <- (start.year + offset):end.year
  } else {
    interval.clim <-(start.year + offset):end.year
    interval.chrono <- (start.year + 1 + offset):end.year
  }
  a <- as.numeric(rownames(chrono)) %in% interval.chrono
  b <- clim[, 1] %in% interval.clim
  chrono.trunc <- chrono[a, 1] # truncate data 
  clim.trunc <- clim[b, ]
  p <- pmat(clim.trunc, start, end, vnames) # call pmat to obtain
                                        # variable matrix
  ## number of windows
  years <- as.numeric(rownames(p))
  n.years <- length(years)
  win.num <- (length(chrono.trunc) - win.size) %/% win.offset
  if (win.num < 2) {
    errmsg6 <-
      paste("Less than 2 windows. Consider a timespan greater than ",
            n.years, " or a win.size smaller than ", win.size, ".",
            sep = "")
    stop(errmsg6)
  }
  win.years.string <- character(win.num)
  windows <- 1:win.num

  ## result matrices
  result.matrix.coef <- result.matrix.ci.upper <-
    result.matrix.ci.lower <- result.matrix.significant <-
      matrix(NA, ncol = win.num, nrow = dim(p)[2])

  if (sb) { # initialize status bar (if TRUE)
    mpb <- txtProgressBar(min = 1,  max = win.num, style = 3)
  }

  METHOD <- match.arg(method, c("response", "correlation"))  # match method argument

  if (startlast) {
    for (k in 1:win.num) {
      series.subset.index <- ((n.years - ((k-1)*win.offset)) -
                              (win.size - 1)):(n.years - ((k-1)*win.offset))
      p.win <- p[series.subset.index,]
      chrono.win <- chrono.trunc[series.subset.index]
      if (METHOD == "response") {
        if (boot) {
          dc.res <- brf(chrono.win, p.win, sb = FALSE, vnames =
                        vnames, ci = ci) # call brf to
                                        # calculate bootstrapped
                                        # response function
          dc.coef <- dc.res$coef
          dc.ci.upper <- dc.res$ci.upper
          dc.ci.lower <- dc.res$ci.lower
          dc.sig <- dc.res$significant
        } else {
          dc.res <- nbrf(chrono.win, p.win, vnames = vnames) # call rf for
                                        # non-bootstrapped response
                                        # function
          dc.coef <- dc.res$coef
          dc.ci.upper <- NA
          dc.ci.lower <- NA
          dc.sig <- NA
        }
      }
      if (METHOD == "correlation") {
        if (boot) {
          dc.res <- bcf(chrono.win, p.win, sb = FALSE, vnames =
                         vnames, ci = ci) # call bcf to
                                        # calculate bootstrapped
                                        # correlation function
          dc.coef <- dc.res$coef
          dc.ci.upper <- dc.res$ci.upper
          dc.ci.lower <- dc.res$ci.lower
          dc.sig <- dc.res$significant
        } else {
          dc.res <- nbcf(chrono.win, p.win, vnames = vnames) # call cf for
                                        # non-bootstrapped correlation
                                        # function
          dc.coef <- dc.res$coef
          dc.ci.upper <- NA
          dc.ci.lower <- NA
          dc.sig <- NA
        }
      }
      result.matrix.coef[,k] <- dc.coef
      result.matrix.ci.upper[,k] <- dc.ci.upper
      result.matrix.ci.lower[,k] <- dc.ci.lower
      result.matrix.significant[,k] <- dc.sig
      win.years.string[k] <- paste(years[series.subset.index][1],
                                   years[series.subset.index][win.size],
                                   sep = "-")
      if (sb) # update status bar (if TRUE)
        setTxtProgressBar(mpb, k)
    }
    result.matrix.coef <- result.matrix.coef[,win.num:1]
    result.matrix.ci.upper <- result.matrix.ci.upper[,win.num:1]
    result.matrix.ci.lower <- result.matrix.ci.lower[,win.num:1]
    result.matrix.significant <- result.matrix.significant[,win.num:1]
    win.years.string <- win.years.string[win.num:1]
  } else {
    for (k in 1:win.num) {
      series.subset.index <- (1 + ((k-1) * win.offset)):(1 + ((k-1) *
                                                              win.offset) +
                                                         (win.size - 1))
      p.win <- p[series.subset.index,]
      chrono.win <- chrono.trunc[series.subset.index]
      if (METHOD == "response") {
        if (boot) {
          dc.res <- brf(chrono.win, p.win, sb = FALSE, vnames =
                         vnames, ci = ci) # call brf to
                                        # calculate bootstrapped
                                        # response function
          dc.coef <- dc.res$coef
          dc.ci.upper <- dc.res$ci.upper
          dc.ci.lower <- dc.res$ci.lower
          dc.sig <- dc.res$significant
        } else {
          dc.res <- nbrf(chrono.win, p.win, vnames = vnames) # call rf for
                                        # non-bootstrapped response
                                        # function
          dc.coef <- dc.res$coef
          dc.ci.upper <- NA
          dc.ci.lower <- NA
          dc.sig <- NA
        }
      }
      if (METHOD == "correlation") {
        if (boot) {
          dc.res <- bcf(chrono.win, p.win, sb = FALSE, vnames =
                         vnames, ci = ci) # call bcf to
                                        # calculate bootstrapped
                                        # correlation function
          dc.coef <- dc.res$coef
          dc.ci.upper <- dc.res$ci.upper
          dc.ci.lower <- dc.res$ci.lower
          dc.sig <- dc.res$significant
        } else {
          dc.res <- nbcf(chrono.win, p.win, vnames = vnames) # call cf for
                                        # non-bootstrapped correlation
                                        # function
          dc.coef <- dc.res$coef
          dc.ci.upper <- NA
          dc.ci.lower <- NA
          dc.sig <- NA
        }
      }
      result.matrix.coef <- result.matrix.coef[,win.num:1]
      result.matrix.ci.upper <- result.matrix.ci.upper[,win.num:1]
      result.matrix.ci.lower <- result.matrix.ci.lower[,win.num:1]
      result.matrix.significant <- result.matrix.significant[,win.num:1]
      win.years.string[k] <- paste(years[series.subset.index][1],
                                   years[series.subset.index][win.size],
                                   sep = "-")
      if (sb) # update status bar (if TRUE)
        setTxtProgressBar(mpb, k)
    }
  }
  dc <- list()
  dc$coef <- data.frame(result.matrix.coef)
  colnames(dc$coef) <- win.years.string
  rownames(dc$coef) <- colnames(p.win)
  dc$ci.upper <- data.frame(result.matrix.ci.upper)
  colnames(dc$ci.upper) <- win.years.string
  rownames(dc$ci.upper) <- colnames(p.win)
  dc$ci.lower <- data.frame(result.matrix.ci.lower)
  colnames(dc$ci.lower) <- win.years.string
  rownames(dc$ci.lower) <- colnames(p.win)
  dc$significant <- data.frame(result.matrix.significant)
  colnames(dc$significant) <- win.years.string
  rownames(dc$significant) <- colnames(p.win)
  if (sb) # close status bar (if TRUE)
    close(mpb)
  dc
}

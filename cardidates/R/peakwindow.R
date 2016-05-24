`peakwindow` <-
function(x, y = NULL, xstart = 0, xmax = max(x),
                        minpeak = 0.1, mincut = 0.382) {
  ## local function for identifying start and end of all peaks
  ## does also numbering of peaks (id-number)
  numpeaks <- function(ft, y) {
    npeaks <- length(ft)
    ndata  <- length(y)
    if (npeaks > 0) {
      fp <- x1 <- x2 <- numeric(npeaks)
      peakid <- numeric(length(y))
      for (i in 1:npeaks) {
        fp[i] <- ft[i]
        mp <- match(fp[i], ft)
        if (mp == 1)      x1[i] <- 1     else x1[i] <- ft[mp-1] + which.min(y[ft[mp-1]:fp[i]])-1
        if (mp == npeaks) x2[i] <- ndata else x2[i] <- fp[i]    + which.min(y[fp[i]:ft[mp+1]])-1
        peakid[x1[i]:x2[i]] <- i
      }
    } else {
      fp <- ft
      x1 <- x2 <- NULL
      peakid <- rep(0, ndata)
    }

    list(fp=fp, x1=x1, x2=x2, id=peakid)
  }

  ## start of main algorithm
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y

  ## Technical note: xmax = max(x) is a "promise".
  ## It is evaluated *now* with the new x after xy.coords
  iend <- max(c(1, which(x <= xmax)))
  tp   <- turnpoints(y)
  ft   <- tp$tppos
  
  ## default return value in case of no turnpoints
  if (all(is.na(ft))) {
    x1 <- 1
    x2 <- length(y)
    fp <- NULL
    ft <- NULL
  } else {
    ## Determination of all peaks
    ft   <- ft[seq(2 - 1 * tp$firstispeak, length(ft), by = 2)]
    ## Add first point if it is a peak
    if (y[1] > y[2]) ft <- c(1, ft)
    ## Add last point if it is a peak
    if (y[iend] > y[iend - 1]) ft <- c(ft, iend)
    ## Limitation until index iend (e.g. iend = mid of year)
    ft   <- ft[ft <= iend]
    ## Limitation to peaks > minpeak * maxvalue (e.g. 10%)
    ftsmall <- which(y < minpeak * max(y))
    ft      <- ft[!ft %in% ftsmall]
    ## If more than one peak then search ...
    dosearch <- ifelse(length(ft) > 1, TRUE, FALSE)
    while(dosearch) {
      fl  <- length(ft)
      eli <- NULL
      for (n in 1:(fl - 1)) {
        ## Minimum between 2 peaks
        km <- which.min(y[ft[n:(n+1)]])
        ## Heuristics:
        ## if Minimum between 2 peaks > mincut of the smaller neighbours
        ## then eliminate peak
        if (min(y[ft[n:(n+1)]]) * mincut < min(y[ft[n]:ft[n+1]]))
        eli <- c(eli, n + km - 1)
      }
      if (length(eli) == 0) dosearch <- FALSE else ft <- ft[-eli]
      if (length(ft) < 2)   dosearch <- FALSE
    }
    ## Search as long as
    ## - no peak can be eliminated or
    ## - only one peak left over
    ## --> ToDo: shorten these lines so that they use "numpeaks()"
    if (length(ft) > 1) {
      ## If more than one peak
      ## take the first one after xstart (day after ice out)
      ## If no peak after ice out than take the last one for xstart
      if (any(x[ft] > xstart)) {
        fp <- min(ft[x[ft]  > xstart])
      } else {
        fp <- max(ft[x[ft] <= xstart])
      }
      mp <- match(fp, ft)
      if (mp == 1)          x1 <- 1  else x1 <- ft[mp-1] + which.min(y[ft[mp-1]:fp])-1
      if (mp == length(ft)) x2 <- iend else x2 <- fp       + which.min(y[fp:ft[mp+1]])-1
    } else {
      fp <- ft
      x1 <- 1
      x2 <- iend
    }
  }
  allpeaks <- numpeaks(ft, y)
  ## Output: - all peak maxima as data frame
  ##         - Start of data window, peak, end of window
  ##           given as indices, x and y values

  ret <-  list(peaks = data.frame(
                         index = ft,
                         xleft = allpeaks$x1,
                         x     = x[ft],
                         xright= allpeaks$x2,
                         y     = y[ft]
                       ),
               data   = data.frame(x=x, y=y),
               smd.max.index = fp,
               smd.max.x     = x[fp],
               smd.indices   = x1:x2,
               smd.x         = x[x1:x2],
               smd.y         = y[x1:x2],
               peakid = allpeaks$id
  )
  class(ret) <- c("list", "cardiPeakwindow")
  ret
}


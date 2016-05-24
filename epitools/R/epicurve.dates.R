"epicurve.dates" <-
  function(x, format = "%Y-%m-%d", strata = NULL,
           min.date, max.date, before = 7, after = 7,
           width = 1, space = 0, tick = TRUE,
           tick.offset = 0.5, segments = FALSE, ...){
    dates0 <- as.Date(x, format = format)
    if(missing(min.date)){
      min.date <- min(dates0, na.rm=TRUE) - before
    }
    if(missing(max.date)){
      max.date <- max(dates0, na.rm=TRUE) + after
    }
    cdates <- seq(min.date, max.date, by = 1)
    dates <- factor(dates0, levels = cdates)
    if(is.null(strata)){
      dat <- t(as.matrix(table(dates)))
    } else {
      dat <- t(table(dates, strata))
    }
    xvals <- barplot(dat, width=width, space=space, ...)
    if(tick){
      axis(1, at=c(0, xvals + tick.offset), labels=FALSE, tick=TRUE)
    }
    if(segments){
      x <- xvals-(width/2)
      y2 <- apply(dat,2,sum)
      xy2 <- cbind(x,y2)
      y0 <- cbind(xy2[1,1],0:xy2[1,2])
      z0 <- cbind(y0, y0[,1]+width, y0[,2])
      for(i in 2:nrow(xy2)){
        yy <- cbind(xy2[i,1],0:xy2[i,2])
        z <- cbind(yy, yy[,1]+width, yy[,2])
        z2 <- rbind(z0,z)
        z0 <- z2
      }
      segments(z0[,1],z0[,2],z0[,3],z0[,4])
    }
    cmday <- as.numeric(format(cdates, format = "%d"))
    cmonth <- format(cdates, format = "%b")
    cyear <- format(cdates, format = "%Y")
    rr <- list(dates = dates0,
               dates2 = dates,
               xvals = xvals,
               cdates = cdates,
               cmday = cmday,
               cmonth = cmonth,
               cyear = cyear
               )
    invisible(rr)
}



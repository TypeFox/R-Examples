"epicurve.weeks" <-
  function(x, format = "%Y-%m-%d", strata = NULL,
           min.date, max.date, before = 7, after = 7,
           width = 1, space = 0, tick = TRUE,
           tick.offset = 0.5, segments = FALSE,
           origin = as.Date("1970-01-01"), sunday = TRUE,
           ...){
    aw <- as.week(x, format = format, 
                  min.date = min.date, max.date = max.date,
                  before = before, after = after,
                  sunday = sunday, origin = origin)
    if(sunday) {firstday <- "Sunday"} else {firstday <- "Monday"}
    original.dates <- aw$dates
    cdates <- aw$cstratum
    dates <- aw$stratum2
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
    rr <- list(dates = original.dates,
               firstday= firstday,
               week = aw$week,
               stratum = aw$stratum,
               stratum2 = aw$stratum2,
               stratum3 = aw$stratum3,
               xvals = xvals,               
               cweek = aw$cweek,
               cstratum = aw$cstratum,
               cstratum2 = aw$cstratum2,
               cmday = aw$cmday,
               cmonth = aw$cmonth,
               cyear = aw$cyear
               )
    invisible(rr)
}


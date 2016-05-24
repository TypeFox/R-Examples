"epicurve.months" <-
  function(x, format = "%Y-%m-%d", strata = NULL,
           min.date, max.date, before = 31, after = 31,
           width = 1, space = 0, tick = TRUE,
           tick.offset = 0.5, segments = FALSE,
           origin = as.Date("1970-01-01"),
           ...){
    am <- as.month(x, format = format, 
                  min.date = min.date, max.date = max.date,
                  before = before, after = after,
                  origin = origin)
    original.dates <- am$dates
    cdates <- am$cstratum
    dates <- am$stratum2
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
               mon = am$mon,
               month = am$month,
               stratum = am$stratum,
               stratum2 = am$stratum2,
               stratum3 = am$stratum3,
               xvals = xvals,
               cmon = am$cmon,
               cmonth = am$cmonth,
               cstratum = am$cstratum,
               cstratum2 = am$cstratum2,
               cmday = am$cmday,
               cyear = am$cyear
               )
    invisible(rr)
}


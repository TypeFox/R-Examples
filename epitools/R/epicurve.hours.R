"epicurve.hours" <-
  function(x, mindt, maxdt,
           strata = NULL, half.hour = FALSE,
           width = 1, space = 0, tick = TRUE,
           tick.offset = 0.5, segments = FALSE,
           ...){
    ah <- as.hour(x, mindt = mindt, maxdt = maxdt,
                  half.hour = half.hour)
    xfactor <- ah$stratum3
    if(is.null(strata)){
      dat <- t(as.matrix(table(xfactor)))
    } else {
      dat <- t(table(xfactor, strata))
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
    rr <- list(ct = ah$ct,
               sec = ah$sec,
               min = ah$min,
               hour = ah$hour,
               hour12 = ah$hour12,
               stratum = ah$stratum,
               stratum2 = ah$stratum2,
               stratum3 = ah$stratum3,
               xvals = xvals,
               cstratum = ah$cstratum,
               cstratum2 = ah$cstratum2,
               csec = ah$csec,
               cmin = ah$cmin,
               chour = ah$chour,
               chour12 = ah$chour12,
               campm = ah$campm,
               campm2 = ah$campm2,
               cweekday = ah$cweekday,
               cwkday = ah$cwkday,
               cmday = ah$cmday,
               cmonth = ah$cmonth,
               cmon = ah$cmon,
               cyear = ah$cyear,
               half.hour = ah$half.hour
               )
    invisible(rr)
}

"epicurve.hours" <-
  function(x, mindt, maxdt,
           strata = NULL, half.hour = FALSE,
           width = 1, space = 0, tick = TRUE,
           tick.offset = 0.5, segments = FALSE,
           ...){
    ah <- as.hour(x, mindt = mindt, maxdt = maxdt,
                  half.hour = half.hour)
    xfactor <- ah$stratum3
    if(is.null(strata)){
      dat <- t(as.matrix(table(xfactor)))
    } else {
      dat <- t(table(xfactor, strata))
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
    rr <- list(ct = ah$ct,
               sec = ah$sec,,
               min = ah$min,
               hour = ah$hour,
               hour12 = ah$hour12,
               stratum = ah$stratum,
               stratum2 = ah$stratum2,
               stratum3 = ah$stratum3,
               xvals = xvals,
               cstratum = ah$cstratum,
               cstratum2 = ah$cstratum2,
               csec = ah$csec,
               cmin = ah$cmin,
               chour = ah$chour,
               chour12 = ah$chour12,
               campm = ah$campm,
               campm2 = ah$campm2,
               cweekday = ah$cweekday,
               cwkday = ah$cwkday,
               cmday = ah$cmday,
               cmonth = ah$cmonth,
               cmon = ah$cmon,
               cyear = ah$cyear,
               half.hour = ah$half.hour
               )
    invisible(rr)
}
"epicurve.hours" <-
  function(x, mindt, maxdt,
           strata = NULL, half.hour = FALSE,
           width = 1, space = 0, tick = TRUE,
           tick.offset = 0.5, segments = FALSE,
           ...){
    ah <- as.hour(x, mindt = mindt, maxdt = maxdt,
                  half.hour = half.hour)
    xfactor <- ah$stratum3
    if(is.null(strata)){
      dat <- t(as.matrix(table(xfactor)))
    } else {
      dat <- t(table(xfactor, strata))
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
    rr <- list(ct = ah$ct,
               sec = ah$sec,
               min = ah$min,
               hour = ah$hour,
               hour12 = ah$hour12,
               stratum = ah$stratum,
               stratum2 = ah$stratum2,
               stratum3 = ah$stratum3,
               xvals = xvals,
               cstratum = ah$cstratum,
               cstratum2 = ah$cstratum2,
               csec = ah$csec,
               cmin = ah$cmin,
               chour = ah$chour,
               chour12 = ah$chour12,
               campm = ah$campm,
               campm2 = ah$campm2,
               cweekday = ah$cweekday,
               cwkday = ah$cwkday,
               cmday = ah$cmday,
               cmonth = ah$cmonth,
               cmon = ah$cmon,
               cyear = ah$cyear,
               half.hour = ah$half.hour
               )
    invisible(rr)
}

inferTi <- function(dateTimes){
  ## make an educated guess as to the frequency of dateTimes and return a ti
  ## object of the same length
  naSpots <- is.na(dateTimes)
  hasNAs <- any(naSpots)
  if(all(naSpots)) stop("dateTimes are all NA")
  dt <- as.POSIXct(as.POSIXlt(dateTimes[!naSpots]))
  dtJul <- floor(jul(dt))
  diffSeconds <- median(diff(unique(sort(unclass(dt)))))
  freq <- round((365.25 * 60*60*24)/diffSeconds)
  
  if(freq > 365){  ## maybe intraday
    if((diffSeconds %% 3600) == 0)  tif <- hourly(diffSeconds / 3600)
    else {
      if((diffSeconds %% 60) == 0)  tif <- minutely(diffSeconds / 60)
      else                          tif <- secondly(diffSeconds)
    }
  }
  else {
    if(freq == 365){
      if(all(between(dayOfWeek(dtJul), 2, 6))) tif <- "business"
      else                                     tif <- "daily"
    }
    else tif <- freq2tif(freq)
  }

  lt <- as.POSIXlt(dt)
  if(isIntradayTif(tif) || (sum(lt$sec + lt$min + lt$hour) == 0))
    dtTi <- ti(dt, tif = tif)
  else
    dtTi <- ti(dt - diffSeconds/2, tif = tif)

  if(freq < 365){
    if(median(abs(jul(dtTi) - dtJul)) > 0.5){ ## could be wrong tif
      maxJul <- max(dtJul) ## the most recent date is most likely correct
      newTif <-
        switch(as.character(freq),
               "52" = tif("wsunday")    + dayOfWeek(maxJul) - 1,
               "26" = tif("bw1sunday")  + dayOfPeriod(maxJul, "bw1sunday") - 1,
               "12" = tif("monthly"),
                "6" = tif("bmdecember") - (month(maxJul) %% 2),
                "4" = tif("qoctober")   + ((2 + month(maxJul)) %% 4), 
                "2" = tif("sannjuly")   + ((5 + month(maxJul)) %% 6),
                "1" = tif("annjanuary") - 1 + month(maxJul),
                NULL)
      if(is.null(newTif))
        stop(paste("Could not infer tif from apparent frequency:", freq))
      else 
        dtTi <- ti(dtJul - (diffSeconds/(2*24*60*60)), tif = newTif)

      if(sum(as.numeric(jul(firstDayOf(dtTi + 1)) == dtJul)) > (length(dtJul)/2))
        dtTi <- dtTi + 1
    }
  }
  ans <- numeric(length(naSpots)) + NA
  ans[!naSpots] <- dtTi
  asTi(ans)
}

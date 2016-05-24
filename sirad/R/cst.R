cst <- 
  function (RefRad, days, lat, extraT=NULL, perce = 3, sepYear = FALSE, stat="median") 
  {
    i <- dayOfYear(days)
    if (is.null(extraT)) extraT <- extrat(i, lat)$ExtraTerrestrialSolarRadiationDaily  
    if (sepYear == TRUE) {
      years <- unique(format(as.Date(days, origin = "1970-01-01"), 
                             "%Y"))
      mtalsALL <- vector()
      for (y in years) {
        wh <- which(format(as.Date(days, origin = "1970-01-01"), 
                           "%Y") == y)
        tals <- RefRad[wh]/extraT[wh]
        stals <- sort(tals)
        ll <- length(stals)
        mtals <- stals[(ll - round(ll * (perce/100))):ll]
        mtalsALL <- c(mtalsALL, mtals)
      }
    }
    if (sepYear == FALSE) {
      tals <- RefRad/extraT
      ll <- length(sort(tals))
      mtalsALL <- sort(tals)[(ll - round(ll * (perce/100))):ll]
      
    }
    if (is.numeric(stat)) rval <- quantile(mtalsALL,probs=stat)
    if (stat=="median") rval <- median(mtalsALL)
    if (stat=="max") rval <- max(mtalsALL)
    rval
  }

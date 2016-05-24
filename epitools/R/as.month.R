"as.month" <-
  function(x, format = "%Y-%m-%d",
           min.date, max.date, before = 31, after = 31,
           origin = as.Date("1970-01-01"), abbreviate = TRUE){
    dates <- as.Date(x, format = format)
    posixlt <- as.POSIXlt(dates)
    mday <- posixlt$mday
    mon <- posixlt$mon +1
    if(abbreviate){
      month <- format(dates, format = "%b")
    } else month <- format(dates, format = "%B")
    jul <- julian(dates, origin = origin)
    stratum <- jul
    md <- c(1:14, 16:31); adj <- c(14:1, (-1:(-16)))
    zz <- cbind(md,adj)
    for(i in 1:nrow(zz)){
      stratum[mday==zz[i,1] & !is.na(mday)] <-
        jul[mday==zz[i,1]  & !is.na(mday)] + zz[i,2]
    }
    if(missing(min.date)) {min.date <- min(dates, na.rm=TRUE) - before}
    if(missing(max.date)) {max.date <- max(dates, na.rm=TRUE) + after}
    cdates <- seq(as.Date(min.date), as.Date(max.date), by = 1)
    cposixlt <- as.POSIXlt(cdates)
    cmday <- cposixlt$mday
    cmon <- cposixlt$mon + 1
    if(abbreviate){
      cmonth <- format(cdates, format = "%b")
    } else cmonth <- format(cdates, format = "%B")    
    cjul <- julian(cdates, origin = origin)
    cstratum <- cjul
    for(i in 1:nrow(zz)){
      cstratum[cmday==zz[i,1] & !is.na(cmday)] <-
        cjul[cmday==zz[i,1]  & !is.na(cmday)] + zz[i,2]
    }    
    repeated <- function(x){
      lx <- length(x)
      y <- rep(FALSE, lx)
      for (i in 2:lx){
        if(x[i]==x[i-1]){
          y[i] <- TRUE
        }
      }
      y
    }
    cmon <- cmon[!repeated(cmon)]
    cmonth <- cmonth[!repeated(cmonth)]
    cstratum <- cstratum[!repeated(cstratum)]
    attr(cstratum, "origin") <- origin
    stratum2 <- factor(stratum, levels = cstratum)
    julian2date <- function(x){
      orig <- as.Date(attributes(x)[[1]])
      jorig <- as.numeric(orig)
      seqdates <- seq(from=orig,to=orig+max(x, na.rm=TRUE),by=1)
      seqjulian <- seq(from=jorig,to=jorig+max(x, na.rm=TRUE),by=1)
      seqdates[x+1]
    }
    stratum3 <- julian2date(stratum)
    cstratum2 <- julian2date(cstratum)
    cmday <- as.numeric(format(cstratum2, format = "%d"))
    cyear <- format(cstratum2, format = "%Y")        
    list(dates = unname(dates),
         mon = unname(mon),
         month = unname(month),
         stratum = stratum,
         stratum2 = stratum2,
         stratum3 = stratum3,
         cmon = unname(cmon),
         cmonth = unname(cmonth),
         cstratum = cstratum,
         cstratum2 = cstratum2,
         cmday = cmday,
         cyear = cyear         
         )
  }


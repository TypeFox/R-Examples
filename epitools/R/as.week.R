"as.week" <-
  function(x, format = "%Y-%m-%d", 
           min.date, max.date, before = 7, after = 7,
           origin = as.Date("1970-01-01"), sunday = TRUE){
    if(sunday) {
      firstday <- "Sunday"
      fday <- "%U"
    } else {
      firstday <- "Monday"
      fday <- "%W"
    }
    dates <- as.Date(x, format = format)
    names(dates) <- as.character(dates)
    wday <- as.POSIXlt(dates)$wday
    jul <- julian(dates, origin = origin)
    week <- format(dates, format = fday)
    stratum <- jul
    if(firstday=="Sunday"){
      stratum[wday==0 & !is.na(wday)] <- jul[wday==0 & !is.na(wday)]+3
      stratum[wday==1 & !is.na(wday)] <- jul[wday==1 & !is.na(wday)]+2
      stratum[wday==2 & !is.na(wday)] <- jul[wday==2 & !is.na(wday)]+1
      stratum[wday==4 & !is.na(wday)] <- jul[wday==4 & !is.na(wday)]-1
      stratum[wday==5 & !is.na(wday)] <- jul[wday==5 & !is.na(wday)]-2
      stratum[wday==6 & !is.na(wday)] <- jul[wday==6 & !is.na(wday)]-3
    }
    if(firstday=="Monday"){
      stratum[wday==1 & !is.na(wday)] <- jul[wday==1 & !is.na(wday)]+3
      stratum[wday==2 & !is.na(wday)] <- jul[wday==2 & !is.na(wday)]+2
      stratum[wday==3 & !is.na(wday)] <- jul[wday==3 & !is.na(wday)]+1
      stratum[wday==5 & !is.na(wday)] <- jul[wday==5 & !is.na(wday)]-1
      stratum[wday==6 & !is.na(wday)] <- jul[wday==6 & !is.na(wday)]-2
      stratum[wday==0 & !is.na(wday)] <- jul[wday==0 & !is.na(wday)]-3
    }    
    if(missing(min.date)) {min.date <- min(dates, na.rm=TRUE) - before}
    if(missing(max.date)) {max.date <- max(dates, na.rm=TRUE) + after}
    mindate <- as.Date(min.date)
    mdtest <- c(as.POSIXlt(mindate)$mon[1], as.POSIXlt(mindate)$mday[1])
    if(mdtest[1]==0 && mdtest[2]<7){
      mindate <- mindate - 7
    }
    cdates <- seq(mindate, as.Date(max.date), by = 1)
    names(cdates) <- as.character(cdates)
    cweek <- cweek2 <- format(cdates, format = fday)
    names(cweek2) <- names(cdates)
    for(i in 1:length(cweek)){
      if(cweek2[i]=="00"){cweek2[i] <- cweek2[i-1]}
    }
    if(mdtest[1]==0 && mdtest[2]<7){
      cdates <- cdates[-c(1:7)]
      cweek <- cweek[-c(1:7)]
      cweek2 <- cweek2[-c(1:7)]
    }    
    week2 <- cweek2[names(dates)]
    cwday <- as.POSIXlt(cdates)$wday
    cjul <- julian(cdates, origin = origin)
    cstratum <- cjul
    if(firstday=="Sunday"){
      cstratum[cwday==0 & !is.na(cwday)] <- cjul[cwday==0 & !is.na(cwday)]+3
      cstratum[cwday==1 & !is.na(cwday)] <- cjul[cwday==1 & !is.na(cwday)]+2
      cstratum[cwday==2 & !is.na(cwday)] <- cjul[cwday==2 & !is.na(cwday)]+1
      cstratum[cwday==4 & !is.na(cwday)] <- cjul[cwday==4 & !is.na(cwday)]-1
      cstratum[cwday==5 & !is.na(cwday)] <- cjul[cwday==5 & !is.na(cwday)]-2
      cstratum[cwday==6 & !is.na(cwday)] <- cjul[cwday==6 & !is.na(cwday)]-3
    }
    if(firstday=="Monday"){
      cstratum[cwday==1 & !is.na(cwday)] <- cjul[cwday==1 & !is.na(cwday)]+3
      cstratum[cwday==2 & !is.na(cwday)] <- cjul[cwday==2 & !is.na(cwday)]+2
      cstratum[cwday==3 & !is.na(cwday)] <- cjul[cwday==3 & !is.na(cwday)]+1
      cstratum[cwday==5 & !is.na(cwday)] <- cjul[cwday==5 & !is.na(cwday)]-1
      cstratum[cwday==6 & !is.na(cwday)] <- cjul[cwday==6 & !is.na(cwday)]-2
      cstratum[cwday==0 & !is.na(cwday)] <- cjul[cwday==0 & !is.na(cwday)]-3
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
    cweek <- cweek2[!repeated(cweek2)]
    cstratum <- cstratum[!repeated(cstratum)]
    attr(cstratum, "origin") <- origin
    stratum <- unname(stratum)
    cstratum <- unname(cstratum)
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
    cmonth <- format(cstratum2, format = "%b")
    cyear <- format(cstratum2, format = "%Y")    
    list(dates = unname(dates),
         firstday = firstday,
         week = unname(week2),         
         stratum = stratum,
         stratum2 = stratum2,
         stratum3 = stratum3,
         cweek = unname(cweek),
         cstratum = cstratum,
         cstratum2 = cstratum2,
         cmday = cmday,
         cmonth = cmonth,
         cyear = cyear
         )
  }

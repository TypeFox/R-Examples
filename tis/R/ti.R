## An not-exported environment used to hold list of frequencies
.tiEnv <- new.env()

## Date representations:
## A ymd is always an integer: year*10000 + month *100 + day
## For jul, c(julianDay, hour, minute, second) is encoded as
##    day + (3600*hour + 60*minute + second)/86400
## For ti (Time Index), c(tif, period) is encoded as tif*1e10 + period
##
## TIFs (Time Index Frequency)
## A tif can be a number between 1001 and 4900, or a string giving the
## frequency name.  tifList() returns a list of the non-intraday frequencies.
## Non-intraday frequency numbers are in the 1001 - 1050 range.
## Intraday frequencies like hourly(3) (every third hour) or secondly(15)
## (every 15'th second) can also be used.  They are encoded as follows:
##      hourly(n)   -->  2000 + n
##      minutely(n) -->  3000 + n
##      secondly(n) -->  4000 + n


asTi  <- function(x) structure(x, class = "ti")
asJul <- function(x) structure(x, class = "jul")

as.jul <- function(x, ...) jul(x, ...)
as.ti  <- function(x, ...) ti(x, ...)

year   <- function(x, ...) ymd(x, ...) %/% 10000
month  <- function(x, ...) (ymd(x, ...) %/% 100) %% 100
day    <- function(x, ...) floor(ymd(x, ...) %% 100)
quarter <- function(x, ...) (month(x, ...) + 2) %/% 3

##From julian date to others
jul2ti <- function(jul, tif){
  tifLen <- length(tif)
  if(tifLen > 1 && length(uniq <- unique(tif)) > 1){
    n <- max(length(jul), tifLen)
    juln <- rep(jul, length.out = n)
    tifn <- rep(tif, length.out = n)
    ans <- asTi(unclass(juln))
    for(u in uniq){
      index <- tifn == u
      ans[index] <- julToTi(juln[index], u)
    }
    return(ans)
  }
  else return(julToTi(jul, tif[1]))
}

jul2ymd <- function(jul){
  jul <- unclass(jul)
  rjul <- jul %/% 1
  j <- as.vector(rjul) - 1721119
  y <- (4*j - 1) %/% 146097
  j <- 4*j - 1 - 146097*y
  d <- j %/% 4
  j <- (4*d + 3) %/% 1461
  d <- 4*d + 3 - 1461*j
  d <- (d + 4) %/% 4
  m <- (5*d - 3) %/% 153
  d <- 5*d - 3 - 153*m
  d <- (d + 5) %/% 5
  y <- 100*y + j
  y <- y + as.numeric(m > 9)
  m <- m + ifelse(m < 10, 3, -9)
  ans <- 10000*y + 100*m + d 
  attributes(ans) <- attributes(jul)
  return(ans)
}

## From ti to others
ti2jul <- function(ti){
  tif <- tif(ti)
  if(length(uniq <- unique(tif)) > 1 || any(is.na(ti))){
    ans <- asJul(unclass(ti))
    for(u in uniq){
      if(!is.na(u)){
        index <- tif == u & !is.na(tif)
        ans[index] <- tiToJul(ti[index])
      }
    }
    return(ans)
  }
  else return(tiToJul(ti))
}

ti2ymd <- function(ti){
  tif <- tif(ti)
  if(length(uniq <- unique(tif)) > 1 || any(is.na(ti))){
    ans <- unclass(ti)
    for(u in uniq){
      if(!is.na(u)){
        index <- tif == u & !is.na(tif)
        ans[index] <- tiToYmd(ti[index])
      }
    }
    return(ans)
  }
  else return(tiToYmd(ti))
}

##From ymd to others
ymd2jul <- function(ymd){
  ucymd <- unclass(ymd)
  seconds <- round(86400 * (ucymd %% 1))
  rYmd <- floor(ucymd + 0.5/86400)
  y <- rYmd %/% 10000
  m <- (rYmd - y*10000) %/% 100
  d <- rYmd %% 100
  y  <- y + (y < 0)
  jy <- y - (m <= 2)
  jm <- m + 1 + 12*(m <= 2)
  jul <- floor(365.25*jy) + floor(30.6001*jm) + d + 1720995;
  addfac <- (2 - floor(0.01*jy) + floor(0.0025*jy))
  i <- (( + 31*(m + 12*y)) >= 588829) & (!is.na(ucymd))
  jul[i] <- jul[i] + addfac[i]
  return(jul + seconds/86400)
}

ymd2ti <- function(ymd, tif){
  tifLen <- length(tif)
  if(length(uniq <- unique(tif)) > 1){
    n <- max(length(ymd), tifLen)
    ymdn <- rep(ymd, length.out = n)
    tifn <- rep(tif, length.out = n)
    ans <- asTi(unclass(ymdn))
    for(u in uniq){
      index <- tifn == u
      ans[index] <- ymdToTi(ymdn[index], u)
    }
    return(ans)
  }
  else return(ymdToTi(ymd, tif[1]))
}

## Object oriented wrappers
## constructors for class 'jul'
is.jul <- function(x) inherits(x, "jul")

jul <- function(x, ...) UseMethod("jul")

jul.jul <- function(x, ...) x

jul.POSIXct <- function(x, ...) jul(POSIXlt(x), ...)

jul.POSIXlt <- function(x, ...){
  j <- (jul(10000*(x$year + 1900) + 100*(x$mon + 1) + x$mday) +
        (x$hour/24) + (x$min/1440) + (x$sec/86400))
  class(j) <- "jul"
  j
}

jul.ssDate <- function(x, ...) jul(18991230) + unclass(x)

jul.ti <- function(x, offset = 1, ...){
  z <- stripClass(x, "ti")
  rx <- x
  oneSecond <- 1/86400
  j1 <- unclass(ti2jul(rx))
  if(!between(offset, 0, 1)) stop("offset must be in [0,1]")
  if(offset == 1)
    z[] <- j1
  else {
    j0 <- unclass(ti2jul(rx - 1))
    if(offset < oneSecond)
      z[] <- j0 + oneSecond
    else{
      z[] <- round((offset*j1 + (1 - offset)*j0)*86400)/86400
    }
  }
  class(z) <- c("jul", oldClass(z))
  return(z)
}

jul.Date <- function(x, ...)
  asJul(unclass(as.vector(x + 2440588)))

jul.IDate <- function(x, ...)
  asJul(unclass(as.vector(x + 2440588)))

jul.default <- function(x, ...){
  if(missing(x))       return(jul.Date(Sys.Date()))
  if(is.character(x))  return(jul.Date(as.Date(x, ...)))
  if(couldBeTi(x))     return(jul(asTi(x), ...))
  if(is.ymd(x))        return(asJul(ymd2jul(x)))
  if(is.time(x))       return(asJul(time2jul(x)))
  else                 return(jul.Date(as.Date(x, ...)))
}

as.character.jul <- function (x, ...) format(x, ...)

as.list.jul <- function(x, ...){
  asClassyList(x, ...)
}

c.jul <- function(..., recursive = F)
  structure(c(unlist(lapply(list(...), unclass))), class = "jul")

diff.jul <- function(x, ...) diff(stripClass(x, "jul"), ...)

min.jul <- function(..., na.rm = F)
  structure(min(unlist(lapply(list(...), unclass)), na.rm = na.rm), class = "jul")

max.jul <- function(..., na.rm = F)
  structure(max(unlist(lapply(list(...), unclass)), na.rm = na.rm), class = "jul")

format.jul <- function(x, ...) format(POSIXlt(x), ...)

print.jul <- function(x, ...){
  ymds <- as.character(floor(ymd(x)))
  hmsList <- hms(x)
  if(sum(unlist(hmsList), na.rm = T) > 0){
    ymds <- paste(ymds,
                  substr(format(100 + hmsList$hour), 2, 3),
                  substr(format(100 + hmsList$min), 2, 3),
                  substr(format(100 + hmsList$sec), 2, 3),
                  sep = ":")
    if(any(naSpots <- is.na(x))) ymds[naSpots] <- as.character(NA)
  }
  names(ymds) <- names(x)
  print(ymds, quote = F, ...)
  cat("class: jul\n")
}

rep.jul <- function(x, times, ...) asJul(NextMethod())

seq.jul <- function(...) asJul(NextMethod())

xtfrm.jul <- function(x) as.numeric(x)

"[.jul" <- function(x, ...) asJul(NextMethod())

Ops.jul <- function (e1, e2){
  if(nargs() == 1){
    if(.Generic == "+") return(e1)
    else stop("operation not defined for jul")
  }
  isJul1 <- is.jul(e1)
  isJul2 <- is.jul(e2)
  if(isJul1 && isJul2){
    validOp <- switch(.Generic, "-" =, "<" = , ">" =, "==" =,
                      "!=" =, "<=" =, ">=" = TRUE, FALSE)
    if(validOp)
      return(do.call(.Generic, list(unclass(e1), unclass(e2))))
    else stop("operation not defined for two juls")
  }
  else {
    if(isJul1) z <- do.call(.Generic, list(unclass(e1), e2))
    else       z <- do.call(.Generic, list(e1, unclass(e2)))
    if(is.numeric(z)) class(z) <- "jul"
    return(z)
  }
}

## 'constructors' for ymd's.  Actually, since there isn't a ymd class, it
## doesn't have constructors.  But you can act as if it did. Or something....
ymd         <- function(x, ...) UseMethod("ymd")
ymd.jul     <- function(x, ...) jul2ymd(x)
ymd.ssDate  <- function(x, ...) ymd(jul(x))

ymd.ti      <- function(x, offset = 1, ...){
  if(!between(offset, 0, 1)) stop("offset must be in [0,1]")
  if(offset == 1)  ti2ymd(x)
  else             ymd(jul(x, offset))
}

ymd.default <- function(x, ...){
  if(missing(x)) jul <- jul()
  else           jul <- jul(x, ...)
  return(jul2ymd(jul))
}

ymdList <- function(x){
  z <- ymd(x)
  list( year = z %/% 10000,
       month = (z %% 10000) %/% 100,
         day = floor(z %% 100))
}

hms <- function(x){
  seconds <- round((unclass(jul(x)) %% 1)*86400)
  hour    <- seconds %/% 3600
  seconds <- seconds %% 3600
  min     <- seconds %/% 60
  sec     <- seconds %% 60
  list(hour = hour, min = min, sec = sec)
}

ymdhms <- function(x) c(ymdList(x), hms(x))

## ssDate (spreadsheet date) class
is.ssDate <- function(x) inherits(x, "ssDate")
as.ssDate <- function(x) structure(x, class = "ssDate")
ssDate <- function(x, ...){
  if(missing(x)) return(ssDate(jul()))
  as.ssDate(jul(x, ...) - jul(18991230))
}
c.ssDate <- function(..., recursive = F)
  structure(c(unlist(lapply(list(...), unclass))), class = "ssDate")

min.ssDate <- function(..., na.rm = F)
  structure(min(unlist(lapply(list(...), unclass)), na.rm = na.rm), class = "ssDate")

max.ssDate <- function(..., na.rm = F)
  structure(max(unlist(lapply(list(...), unclass)), na.rm = na.rm), class = "ssDate")

print.ssDate <- function(x, ...){
  print(ymd(x), ...)
  cat("class: ssDate\n")
}

rep.ssDate   <- function(x, times, ...) as.ssDate(NextMethod())
seq.ssDate   <- function(...) as.ssDate(NextMethod())
xtfrm.ssDate <- function(x) as.numeric(x)
"[.ssDate"   <- function(x, ...) as.ssDate(NextMethod())

Ops.ssDate <- function(e1, e2){
  if(nargs() == 1){
    if(.Generic == "+") return(e1)
    else stop("operation not defined for ssDate")
  }
  is.ssDate1 <- is.ssDate(e1)
  is.ssDate2 <- is.ssDate(e2)
  if(is.ssDate1 && is.ssDate2){
    validOp <- switch(.Generic, "-" =, "<" = , ">" =, "==" =,
                      "!=" =, "<=" =, ">=" = TRUE, FALSE)
    if(validOp)
      return(do.call(.Generic, list(unclass(e1), unclass(e2))))
    else stop("operation not defined for two ssDates")
  }
  else {
    if(is.ssDate1) z <- do.call(.Generic, list(unclass(e1), e2))
    else           z <- do.call(.Generic, list(e1, unclass(e2)))
    if(is.numeric(z)) class(z) <- "ssDate"
    return(z)
  }
}

## ti class
is.ti <- function(x) inherits(x, "ti")

ti <- function(x, ...) UseMethod("ti")
ti.jul <- function(x, tif = NULL, freq = NULL,
                   hour = 0, minute = 0, second = 0, ...){
  if(is.null(tif)){
    if(is.null(freq))
      stop("'tif' and 'freq' cannot both be NULL if 'x' is not a ti")
    else tif <- freq2tif(freq)
  }
  intraday <- isIntradayTif(tif)
  if(any(intraday)){
    if(!(missing(hour) && missing(minute) && missing(second)))
      x[intraday] <- (floor(x + .5/86400) + (3600*hour + 60*minute + second)/86400)[intraday]
  }
  return(jul2ti(x, tif))
}

ti.POSIXlt <- function(x, tif, ...){
  nTif <- tif(tif)
  j <- jul(10000*(x$year + 1900) + 100*(x$mon + 1) + x$mday)
  intraday <- isIntradayTif(nTif)
  if(any(intraday)){
    if(all(intraday))
      j <- j + (x$hour/24) + (x$min/1440) + (x$sec/86400)
    else
      j[intraday] <- (j + (x$hour/24) + (x$min/1440) + (x$sec/86400))[intraday]
  }
  ti(j, nTif)
}

ti.POSIXct <- function(x, tif, ...) ti(POSIXlt(x), tif, ...)

ti.ssDate <- function(x, ...) ti(jul(x), ...)

ti.ti <- function(x, tif = NULL, freq = NULL, ...){
  if(is.null(tif)){
    if(is.null(freq)) return(x)
    else tif <- freq2tif(freq)
  }
  return(ti(jul(x), tif, ...))
}

ti.Date <- function(x, ...) ti(jul(x), ...)

ti.default <- function(x, tif = NULL, freq = NULL, ...){
  if(is.null(tif)){
    if(is.null(freq))
      stop("'tif' and 'freq' cannot both be NULL if 'x' is not a ti")
    else tif <- freq2tif(freq)
  }
  if(missing(x))              return(ti(jul(Sys.Date()), tif, ...))
  if(is.null(x))              stop("NULL argument to ti function")
  if(is.character(x))         return(ti(as.Date(x, ...), tif = tif))
  if(couldBeTi(x, tif = tif)) return(ti(asTi(x), tif = tif, ...))
  
  if(is.ymd(x)){
    if(isIntradayTif(tif)) return(ti(asJul(ymd2jul(x)), tif, ...))
    else                   return(ymd2ti(round(x), tif))
  }
  if(is.time(x)){
    if(isIntradayTif(tif)) return(ti(asJul(time2jul(x)), tif, ...))
    else                   return(ti(asJul(floor(time2jul(x))), tif))
  }
  if(length(x) == 2 && is.time(x[1]) && between(x[2],1,1e10 - 1)){
    firstTi <- ti(ISOdatetime(year = floor(x[1]), 1, 1, 0, 0, 0), tif = tif)
    return(firstTi + x[2] - 1)
  }
  else return(ti(as.Date(x, ...), tif))
}

couldBeTi <- function(x, tif = NULL){
  perhaps <- is.numeric(x) && (all(is.finite(x)) & between(x, 1e13, 5e13))
  if(!perhaps) return(FALSE)
  if(is.null(tif)) return(TRUE)
  else {
    nTif <- tif(tif)
    return(all(between(x, 1e10*nTif, 1e10*(nTif + 1))))
  }
}

period <- function(z){
  if(is.ti(z) || couldBeTi(z))
    stripClass(z, "ti") %% 1e10
  else NULL
}

basePeriod <- function(x){
  nTif <- tif(x)
  asTi(nTif * 1e10)
}

as.character.ti <- function(x, ...){
  j <- jul(x)
  ymds <- as.character(ymd(floor(j)))
  intraday <- isIntradayTif(tif(x))
  if(any(intraday)){
    hmsList <- hms(j[intraday])
    ymds[intraday] <- paste(ymds[intraday],
                            substr(format(100 + hmsList$hour), 2, 3),
                            substr(format(100 + hmsList$min), 2, 3),
                            substr(format(100 + hmsList$sec), 2, 3),
                            sep = ":")
  }
  names(ymds) <- names(x)
  ymds
}

as.list.ti <- function(x, ...){
  asClassyList(x, ...)
}
c.ti <- function(..., recursive = F)
  structure(c(unlist(lapply(list(...), unclass))), class = "ti")

min.ti <- function(..., na.rm = F)
  structure(min(unlist(lapply(list(...), unclass)), na.rm = na.rm), class = "ti")

max.ti <- function(..., na.rm = F)
  structure(max(unlist(lapply(list(...), unclass)), na.rm = na.rm), class = "ti")

diff.ti <- function(x, ...) diff(stripClass(x, "ti"), ...)

cycle.ti <- function(x, ...){
  firstTi <- ti(ISOdatetime(year(x),1,1,0,0,0), tif = tif(x))
  x - firstTi + 1
}

frequency.ti <- function(x, ...){
  f <- round(200/(time(x+100) - time(x-100)))
  if(any(index <- f > 100)) ## handle leap years 
    f[index] <- round(2000/(time(x[index]+1000) - time(x[index]-1000)))
  return(f)
}

deltat.ti <- function(x, ...) 1/frequency(x)

format.ti <- function(x, ..., tz = ""){
  z <- stripClass(x, "ti")
  intraday <- isIntradayTif(tif(x))
  if(any(intraday))
    z[intraday] <- format(POSIXct(x[intraday]), tz = tz, ...)
  if(any(!intraday))
    z[!intraday] <- format(POSIXlt(x[!intraday]), ...)
  z
}

print.ti <- function(x, class = TRUE, ...){
  print(as.character(x), quote = FALSE, ...)
  if(class) cat("class: ti\n")
}

rep.ti <- function(x, times, ...) asTi(NextMethod())

seq.ti <- function(...) asTi(NextMethod())

xtfrm.ti <- function(x) as.numeric(x)

"[.ti" <- function(x, ...) asTi(NextMethod())

Ops.ti <- function (e1, e2){
  if(nargs() == 1){
    if(.Generic == "+") return(e1)
    else stop("operation not defined for ti")
  }
  isTi1 <- is.ti(e1)
  isTi2 <- is.ti(e2)
  if(isTi1 && isTi2){
    n <- max(length(e1), length(e2))
    tif1 <- rep(tif(e1), length.out = n)
    tif2 <- rep(tif(e2), length.out = n)
    if(any(tif1 != tif2) && .Generic != "==") stop("different tif\'s")
    validOp <- switch(.Generic, "-" =, "<" = , ">" =, "==" =,
                      "!=" =, "<=" =, ">=" = TRUE, FALSE)
    if(validOp)
      return(do.call(.Generic, list(unclass(e1), unclass(e2))))
    else
      stop("operation not defined for two ti\'s")
  }
  else {
    if(isTi1) ti <- do.call(.Generic, list(unclass(e1), e2))
    else      ti <- do.call(.Generic, list(e1, unclass(e2)))
    class(ti) <- "ti"
    return(ti)
  }
}

## compatibility with R date classes
format.POSIXlt <- function (x, format = "", usetz = FALSE, ...){
  ## added formats:
  ##    "%N"    first letter of month name
  ##    "%q"    quarter number
  if(!inherits(x, "POSIXlt"))  stop("wrong class")
  if(format == ""){
    times <- unlist(unclass(x)[1:3])
    format <- if(all(times[!is.na(times)] == 0)) "%Y-%m-%d"
    else "%Y-%m-%d %H:%M:%S"
  }
  format <- gsub("%N", "@N", gsub("%q", "@q", format))
  
  tmp <- base::format.POSIXlt(x, format, usetz)
  
  if(length(grep("@q", format))){
    qtr <- (as.numeric(base::format.POSIXlt(x, "%m", usetz)) + 2) %/% 3
    for(i in 1:length(tmp))
    tmp[i] <- gsub("@q", qtr[i], tmp[i])
  }
  
  if(length(grep("@N", format))){
    month <- as.numeric(base::format.POSIXlt(x, "%m", usetz))
    for(i in 1:length(tmp))
      tmp[i] <- gsub("@N", substring(month.abb[month[i]], 1, 1), tmp[i])
  }
  tmp
}

format.POSIXct <- function (x, format = "", tz = "", usetz = FALSE, ...){
  if(!inherits(x, "POSIXct")) stop("wrong class")
  if(missing(tz) && !is.null(tzone <- attr(x, "tzone"))) tz <- tzone
  structure(format(POSIXlt(x, tz), format, usetz, ...),
            names = names(x))
}

weekdays.default <- function(x, ...) weekdays(as.Date(x), ...)
months.default   <- function(x, ...) months(as.Date(x), ...)
quarters.default <- function(x, ...) quarters(as.Date(x), ...)

as.Date.jul <- function(x, ...) structure(x - 2440588, class = "Date")
as.Date.ti  <- function(x, ...) as.Date(jul(x), ...)

as.POSIXct.jul <- function(x, tz = "", ...){
  do.call("ISOdatetime", c(ymdhms(x), tz = tz))
}

as.POSIXct.ti <- function(x, tz = "", offset = 1, ...){
  if(!between(offset, 0, 1)) stop("offset must be in [0,1]")
  
  tiToPOSIXct <- function(x, hour = 23, min = 59, sec = 59){
    argList <- ymdList(x)
    intra <- isIntradayTif(tif(x))
    hmsList <- list(hour = hour, min = min, sec = sec)
    if(any(intra)){
      hmsList <- hms(x)
      if(!all(intra)){
        notIntra <- !intra
        hmsList$hour[notIntra] <- hour
        hmsList$min[notIntra]  <- min
        hmsList$sec[notIntra]  <- sec
      }
    }
    do.call("ISOdatetime", c(ymdList(x), hmsList, tz = tz))
  }
  
  z <- stripClass(x, "ti")
  ct1 <- unclass(tiToPOSIXct(x))
  if(offset == 1) z[] <- ct1
  else {
    ct0 <- unclass(tiToPOSIXct(x - 1) + 1)
    z[] <- floor(offset*(ct1 + 1) + (1 - offset)*ct0)
  }
  attr(z, "tzone") <- attr(ct1, "tzone")
  class(z) <- c("POSIXt", "POSIXct", oldClass(z))
  return(z)
}

POSIXct         <- function(x, ...) UseMethod("POSIXct")
POSIXct.jul     <- function(x, ...) as.POSIXct(x, ...)
POSIXct.ti      <- function(x, offset = 1, ...) as.POSIXct(x, offset = offset, ...)
POSIXct.POSIXt  <- function(x, ...) as.POSIXct(x, ...)
POSIXct.numeric <- function (x, tz = "", origin, ...){
  if(missing(origin)){
    jd <- NULL
    if(couldBeTi(x)) jd <- jul(asTi(x), ...)
    if(is.ymd(x))    jd <- asJul(ymd2jul(x))
    if(is.time(x))   jd <- asJul(time2jul(x))
    
    if(!is.null(jd)) as.POSIXct(jd, tz = tz)
    else             stop("'origin' must be supplied")
  }
  else as.POSIXct(origin, tz = tz, ...) + x
}
POSIXct.default <- function(x, ...) POSIXct(jul(x), ...)

POSIXlt         <- function(x, ...) UseMethod("POSIXlt")
POSIXlt.jul     <- function(x, ...) as.POSIXlt(POSIXct(x, ...))
POSIXlt.ti      <- function(x, ...) as.POSIXlt(POSIXct(x, ...))
POSIXlt.POSIXt  <- function(x, ...) as.POSIXlt(x, ...)
POSIXlt.default <- function(x, ...) POSIXlt(jul(x), ...)

as.POSIXlt.jul  <- function(x, ...) POSIXlt(x, ...)
as.POSIXlt.ti   <- function(x, ...) POSIXlt(x, ...)

baseYmd <- function(tif){
  nTif <- tif(tif)
  if(between(nTif, 1001, 1050)){
    c(## daily, business day
      18991230, 18991229,
      ## wsunday thru wsaturday
      18991218:18991224,
      ## tenday
      18991221,
      ## reserves
      18991221,
      ## biweekly1Sunday thru biweekly2Saturday
      18991211:18991224,
      ## semimonthly and monthly
      18991215, 17991201,
      ## bimonth1 and bimonth2
      rep(17991101, 2),
      ## quarterly
      rep(17991001, 3),
      ## annjanuary thru annnovember
      100*(159802:159812) + 1,
      ## anndecember
      15990101,
      ## sannjuly thru sannnovember
      rep(15990701, 6))[nTif - 1000]      
  }
  else {
    if(between(nTif, 2000, 4900)) 19800101
    else stop(paste("unknown nTif:", nTif))
  }
}

## workhorse functions that convert ti's back and forth to ymd and jul
julToTi <- function(jul, tif, must.handle=F){
  nTif <- tif(tif)
  j <- unclass(jul)
  if(any(jul2ymd(jul) < baseYmd(nTif)))
    stop("jul date too early for tif")
  
  halfSecond <- 1/(86400*2)
  if(!isIntradayTif(nTif))
    j <- floor(j + halfSecond)

      
  if(between(nTif, 1001, 1009) || between(nTif, 1011, 1025)){
    
    period <- switch(nTif - 1e3,  ## handle day-based freqs
                     ## 1 = daily
                     j - 2415019,
                     ## 2 = business day
                     { 
                       dow <- julToWeekday(j)
                       j <- j + (dow==1) + 2*(dow==7)
                       ((j - 2415021)%/%7)*5 + (j - 2415020)%%7
                     },
                     ## 3 - 9 = weeklySunday thru weeklySaturday
                     (j - 2415007)%/%7,
                     (j - 2415008)%/%7,
                     (j - 2415009)%/%7,
                     (j - 2415010)%/%7,
                     (j - 2415011)%/%7,
                     (j - 2415012)%/%7,
                     (j - 2415013)%/%7,
                     ## 10 is tenday, diverted by the if statement above 
                     NULL,
                     ## 11 (reserves) is weekly Wednesday until 19840201 and
                     ## biweekly Wednesday afterwards. firstCrrDay is 19840202.
                     {
                        firstCrrDay <- 2445733
                        preCrr <- j < firstCrrDay
                        prePeriods  <-  (j - 2415010)%/%7
                        postPeriods <- 4389 + (j - firstCrrDay)%/%14
                        prePeriods*preCrr + postPeriods*(!preCrr)
                     },
                     ## 12 - 25 = biweekly1Sunday thru biweekly2Saturday
                     (j - 2415000)%/%14,
                     (j - 2415001)%/%14,
                     (j - 2415002)%/%14,
                     (j - 2415003)%/%14,
                     (j - 2415004)%/%14,
                     (j - 2415005)%/%14,
                     (j - 2415006)%/%14,
                     (j - 2415007)%/%14,
                     (j - 2415008)%/%14,
                     (j - 2415009)%/%14,
                     (j - 2415010)%/%14,
                     (j - 2415011)%/%14,
                     (j - 2415012)%/%14,
                     (j - 2415013)%/%14)
    ans <- asTi(1e10*nTif + period)
    names(ans) <- names(jul)
    return(ans)
  }
  if(isIntradayTif(nTif)){
    j19800101 <- 2444240    ## julian date for Jan 1, 1980
    if(any(j < j19800101)) stop("Intraday ti cannot be earlier than 19800101")
    hms <- nTif %/% 1e3
    nUnits <- nTif %% 1e3
    if(nUnits == 0) nUnits <- 1
    period <- switch(hms - 1, ## 1 = hourly, 2 = minutely, 3 = secondly
                     1 + floor((j + halfSecond - j19800101)*(24/nUnits)),
                     1 + floor((j + halfSecond - j19800101)*(1440/nUnits)),
                     1 + floor((j + halfSecond - j19800101)*(86400/nUnits)))
    ans <- asTi(1e10*nTif + period)
    names(ans) <- names(jul)
    return(ans)
  }  
  else { ## call ymdToTi for other frequencies
    if(must.handle)
      stop(paste(tif, "is an unknown frequency"))
    else 
      return(do.call("ymdToTi",
                     list(ymd = jul2ymd(j), tif = tif, must.handle=T)))
    
  }
}

ymdToTi <- function(ymd, tif, must.handle=F){
  rawYmd <- unclass(ymd)
  year   <- rawYmd %/% 10000
  month  <- (rawYmd - year*10000) %/% 100
  day    <- rawYmd %% 100
  
  nTif <- tif(tif)
  if(any(rawYmd < baseYmd(nTif)))
    stop("ymd too early for tif")
  if(nTif == 1010){ ## tenday
    period <- 1 + 36*(year-1900) + 3*(month-1) + (day>10) + (day>20)
    ans <- asTi(1e10*nTif + period)
    names(ans) <- names(ymd)
    return(ans)
  }
  if(between(nTif, 1026, 1050)){
    period <- switch(nTif - 1e3 - 25,
                     ## 26 = twicemonthly
                     1 + 24*(year-1900) + 2*(month-1) + (day>15),
                     ## 27 = monthly
                     12*(year-1800) + month,
                     ## 28 and 29 are bimonthly (nov and dec)
                     6*(year-1800) + (month+2) %/% 2, 
                     6*(year-1800) + (month+1) %/% 2,
                     ## 30 - 32 are quarterly (oct, nov and dec)
                     4*(year-1800) + (month+4) %/% 3, 
                     4*(year-1800) + (month+3) %/% 3, 
                     4*(year-1800) + (month+2) %/% 3,
                     ## 33 - 44 are annual (jan, feb, ..., dec)
                     1 + (year-1600) + (month > 1),
                     1 + (year-1600) + (month > 2),
                     1 + (year-1600) + (month > 3),
                     1 + (year-1600) + (month > 4),
                     1 + (year-1600) + (month > 5),
                     1 + (year-1600) + (month > 6),
                     1 + (year-1600) + (month > 7),
                     1 + (year-1600) + (month > 8),
                     1 + (year-1600) + (month > 9),
                     1 + (year-1600) + (month > 10),
                     1 + (year-1600) + (month > 11),
                     1 + (year-1600),
                     ## 45 - 50 are semi-annual (jul, aug, ..., dec)
                     1 + 2*(year-1600) + (month > 6) + ((month %% 6) > 1),
                     1 + 2*(year-1600) + (month > 6) + ((month %% 6) > 2),
                     1 + 2*(year-1600) + (month > 6) + ((month %% 6) > 3),
                     1 + 2*(year-1600) + (month > 6) + ((month %% 6) > 4),
                     1 + 2*(year-1600) + (month > 6) + ((month %% 6) > 5),
                     1 + 2*(year-1600) + (month > 6))
    ans <- asTi(1e10*nTif + period)
    names(ans) <- names(ymd)
    return(ans)
  }
  else { ## call ti.jul for other frequencies
    if(must.handle)
      stop(paste(tif, "is an unknown frequency"))
    else 
      return(do.call("julToTi",
                     list(jul = ymd2jul(ymd),tif = tif, must.handle=T)))
  }
}

tiToJul <- function(ti, must.handle=F){
  uti <- as.vector(unclass(ti))
  nTif <- uti[1] %/% 1e10
  periods <- (uti %% 1e10) - 1
  
  if(between(nTif, 1001, 1009) || between(nTif, 1011, 1025)){
    j <- switch(nTif - 1e3,  ## handle day-based freqs
                ## 1 = daily
                2415020 + periods,
                ## 2 = business day
                2415021 + 7*(periods %/% 5) + periods %% 5,
                ## 3 - 9 = weeklySunday thru weeklySaturday
                2415020 + periods*7,
                2415021 + periods*7,
                2415022 + periods*7,
                2415023 + periods*7,
                2415024 + periods*7,
                2415025 + periods*7,
                2415026 + periods*7,
                ## Use placeholder for 10 (= tenday)
                NULL,
                ## need code for 11 (= reserves)
                {
                  preCrr  <- periods < 4388
                  preJul  <- 2415023 + periods*7
                  postJul <- 2445732 + (periods - 4387)*14
                  preJul*preCrr + postJul*(!preCrr)
                },
                ## 12 - 25 = biweekly1Sunday thru biweekly2Saturday
                2415027 + periods*14,
                2415028 + periods*14,
                2415029 + periods*14,
                2415030 + periods*14,
                2415031 + periods*14,
                2415032 + periods*14,
                2415033 + periods*14,
                2415034 + periods*14,
                2415035 + periods*14,
                2415036 + periods*14,
                2415037 + periods*14,
                2415038 + periods*14,
                2415039 + periods*14,
                2415040 + periods*14)
    names(j) <- names(ti)
    return(asJul(j))
  }
  else{
    if(isIntradayTif(nTif)){
      j19800101 <- 2444240    ## julian date for Jan 1, 1980
      hms <- nTif %/% 1e3
      nUnits <- nTif %% 1e3
      if(nUnits == 0) nUnits <- 1
      j <- switch(hms - 1, ## 1 = hourly, 2 = minutely, 3 = secondly
                  j19800101 + nUnits*periods/24,
                  j19800101 + nUnits*periods/1440,
                  j19800101 + nUnits*periods/86400)
      names(j) <- names(ti)
      return(asJul(j))
    }
    else { ## pass other frequencies to tiToYmd
      if(must.handle)  stop(paste(nTif, "is an unknown frequency"))
      else
        return(ymd2jul(tiToYmd(ti, must.handle=T)))
    }
  }
}

tiToYmd <- function(ti, must.handle=F){
  uti <- as.vector(unclass(ti))
  nTif <- uti[1] %/% 1e10
  periods <- (uti %% 1e10) - 1
  mdays <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  if(nTif == 1010){ ## tenday
    y <- 1900 + periods %/% 36
    m <- 1 + (periods %% 36) %/% 3
    d <- mdays[m]
    third <- periods %% 3
    d[third == 0] <- 10
    d[third == 1] <- 20
    d <- d + (isLeapYear(y) & m == 2 & d > 27)
    ymd <- y*10000 + m*100 + d
    names(ymd) <- names(ti)
    return(ymd)
  }
  else {
    if(between(nTif, 1026, 1050)){
      switch(nTif - 1e3 - 25, 
             ## 26 = twicemonthly
             { y <- 1900 + periods %/% 24
               m <- 1 + (periods %% 24) %/% 2
               half <- periods %% 2
               d <- 15*(half==0) + mdays[m]*(half==1)
             },
             ## 27 = monthly
             { y <- 1800 + periods %/% 12
               m <- 1 + periods %% 12
               d <- mdays[m]
             },
             ## 28 and 29 are bimonthly (nov and dec)
             { y <- 1800 + periods %/% 6
               m <- 1 + 2*(periods %% 6)
               d <- mdays[m]
             }, 
             { y <- 1800 + periods %/% 6
               m <- 2 + 2*(periods %% 6)
               d <- mdays[m]
             },
             ## 30 - 32 are quarterly (oct, nov and dec)
             { y <- 1800 + periods %/% 4
               m <- 1 + 3*(periods %% 4)
               d <- mdays[m]
             },
             { y <- 1800 + periods %/% 4
               m <- 2 + 3*(periods %% 4)
               d <- mdays[m]
             }, 
             { y <- 1800 + periods %/% 4
               m <- 3 + 3*(periods %% 4)
               d <- mdays[m]
             },
             ## 33 - 44 are annual (jan, feb, ..., dec)
             { y <- 1600 + periods
               m <- 1
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 2
               d <- mdays[m]
             },
             { y <- 1600 + periods
               m <- 3
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 4
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 5
               d <- mdays[m]
             },
             { y <- 1600 + periods
               m <- 6
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 7
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 8
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 9
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 10
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 11
               d <- mdays[m]
             }, 
             { y <- 1600 + periods
               m <- 12
               d <- mdays[m]
             },
             ## 45 - 50 are semi-annual (jul, aug, ..., dec)
             { y <- 1600 + periods %/% 2
               m <- 1 + 6*(periods %% 2)
               d <- mdays[m]
             }, 
             { y <- 1600 + periods %/% 2
               m <- 2 + 6*(periods %% 2)
               d <- mdays[m]
             }, 
             { y <- 1600 + periods %/% 2
               m <- 3 + 6*(periods %% 2)
               d <- mdays[m]
             }, 
             { y <- 1600 + periods %/% 2
               m <- 4 + 6*(periods %% 2)
               d <- mdays[m]
             }, 
             { y <- 1600 + periods %/% 2
               m <- 5 + 6*(periods %% 2)
               d <- mdays[m]
             }, 
             { y <- 1600 + periods %/% 2
               m <- 6 + 6*(periods %% 2)
               d <- mdays[m]
             })
    }
    else { ## pass other frequencies to tiToJul
      if(must.handle)
        stop(paste(tif, "is an unknown frequency"))
      else
        return(jul2ymd(tiToJul(ti, must.handle=T)))
    }
  }
  d <- d + (isLeapYear(y) & m == 2 & d > 27)
  ymd <- y*10000 + m*100 + d
  names(ymd) <- names(ti)
  return(ymd)
}

tif <- function(x, ...) UseMethod("tif")
tif.ti  <- function(x, ...) unclass(x) %/% 1e10
tif.tis <- function(x, ...) tif(start(x))
tif.ts  <- function(x, ...) tif(as.tis(x))
tif.default <- function(x, freq = NULL, ...){
  ## freq ignored unless missing(x), ... is always ignored
  if(missing(x)){
    if(is.null(freq)) return(tifList())
    else return(freq2tif(freq))
  }
  if(is.numeric(x) && between(x, 1001, 4900))
    return(as.vector(unclass(x)))
  else {
    if(is.character(x)){
      xlen <- length(x)
      ans <- numeric(xlen) + NA
      intraday <- substr(x, 1, 6) %in% c("hourly", "minute", "second")
      if(any(intraday)){
        spots <- (1:xlen)[intraday]
        for(spot in spots) ans[spot] <- eval(parse(text = x[spot]))
      }
      if(any(!intraday))
        ans[!intraday] <- as.vector(tifList()[x[!intraday]])
      if(any(is.na(ans))) stop("no such tif")
      else return(ans)
    }
    else stop("arg must be ti, tis, ts, tif, or tifName")
  }
}

## Intraday frequencies
isIntradayTif <- function(tif){
  if(is.character(tif)) tif <- tif(tif)
  z <- between(tif, 2000, 4900)
  z[is.na(z)] <- F
  z
}
hourly <- function(n = 0){
  if(n == 0) return(2000)
  if(!between(n, 0, 24)) stop("n out of range")
  if(!(24 %% n == 0)) stop("n not a factor of 24")
  return(2000 + n)
}
minutely <- function(n = 0){
  if(n == 0) return(3000)
  if(!between(n, 0, 1024)) stop("n out of range")
  if((n %% 60) == 0) return(hourly(n %/% 60))
  if(!(1440 %% n == 0)) stop("n not a factor of 1440")
  return(3000 + n)
}
secondly <- function(n = 0){
  if(n == 0) return(4000)
  if(!between(n, 0, 1024)) stop("n out of range")
  if((n %% 60) == 0) return(minutely(n %/% 60))
  if(!(86400 %% n == 0)) stop("n not a factor of 86400")
  return(4000 + n)
}

## Support functions
isLeapYear    <- function(y) y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0)
is.ymd        <- function(x) all(between(x, 15830101, 29991231), na.rm = T)
julToWeekday  <- function(jul){
  ## Sun = 1, Sat = 7.  2415020 = Sunday, 12/31/1899
  ((unclass(jul) - 2415020) %% 7) + 1
}

freq2tif <- function(freq){
  if(is.null(freq)) stop("NULL freq")
  if(!is.numeric(freq)) stop("freq must be numeric")
  tif <- switch((freq >= 1) + (freq >= 2) + (freq >= 4) + 
                (freq >= 6) + (freq >= 12) + (freq >= 24) + 
                (freq >= 26) + (freq >= 36) + (freq >= 52) + 
                (freq >= 261) + (freq >= 365), 
                44, 50, 32, 29, 27, 26, 22, 10, 4, 2, 1)
  if(is.null(tif)) stop("can't convert freq to tif")
  else             return(1e3 + tif)
}

tif2freq <- function(tif)  frequency(ti(tif = tif))

tifName <- function(s) UseMethod("tifName")
tifName.tis <- function(s) tifName(start(s))
tifName.ti <- function(s){
  nTif <- tif(s)
  cTif <- stripClass(s, "ti")
  mode(cTif) <- "character"
  intraday <- isIntradayTif(nTif)
  if(any(intraday)){
    base <- c("hourly", "minutely", "secondly")[(nTif[intraday] %/% 1000) -1]
    nUnits <- as.character(nTif[intraday] %% 1000)
    nUnits[nUnits == "0"] <- ""
    cTif[intraday] <- paste(base, "(", nUnits, ")", sep = "")
  }
  if(any(notIntraday <- !intraday)){
    tl <- tifList()
    cTif[notIntraday] <- names(tl)[match(nTif[notIntraday], tl)]
  }
  cTif
}

tifName.default <- function(s){
  tl <- tifList()
  if(missing(s)) return(tl)
  if(is.character(s)){
    if(all((s %in% names(tl)) |
           (substr(s, 1, 6) %in% c("hourly", "minute", "second"))))
      return(s)
    else stop("unknown tifName")
  }
  else{
    if(!is.numeric(s)) stop("non-numeric, non-character arg s")
    cTif <- character(length(s))
    intraday <- isIntradayTif(s)
    if(any(intraday)){
      base <- c("hourly", "minutely", "secondly")[(s[intraday] %/% 1000) - 1]
      nUnits <- as.character(s[intraday] %% 1000)
      nUnits[nUnits == "0"] <- ""
      cTif[intraday] <- paste(base, "(", nUnits, ")", sep = "")
    }
    if(any(notIntraday <- !intraday))
      cTif[notIntraday] <- names(tl)[match(s[notIntraday], tl)]
    return(cTif)
  }
}

initialTifList <- function(){
  c(daily               = 1001,
    business            = 1002,
    wsunday             = 1003,
    wmonday             = 1004,
    wtuesday            = 1005,
    wwednesday          = 1006,
    wthursday           = 1007,
    wfriday             = 1008,
    wsaturday           = 1009,
    tenday              = 1010,
    reserves            = 1011,
    bw1sunday           = 1012,
    bw1monday           = 1013,
    bw1tuesday          = 1014,
    bw1wednesday        = 1015,
    bw1thursday         = 1016,
    bw1friday           = 1017,
    bw1saturday         = 1018,
    bw2sunday           = 1019,
    bw2monday           = 1020,
    bw2tuesday          = 1021,
    bw2wednesday        = 1022,
    bw2thursday         = 1023,
    bw2friday           = 1024,
    bw2saturday         = 1025,
    twicemonth          = 1026,
    semimonthly         = 1026,
    monthly             = 1027,
    bimonth1            = 1028,
    bmnovember          = 1028,
    bimonthnovember     = 1028,
    bimonth2            = 1029,
    bmdecember          = 1029,
    bimonthdecember     = 1029, 
    bimonth             = 1029,
    quarterlyoctoctober = 1030,
    qoctober            = 1030,
    quarterlynovember   = 1031,
    qnovember           = 1031,
    quarterlydecember   = 1032,
    qdecember           = 1032,
    annjanuary          = 1033,
    annfebruary         = 1034,
    annmarch            = 1035,
    annapril            = 1036,
    annmay              = 1037,
    annjune             = 1038,
    annjuly             = 1039,
    annaugust           = 1040,
    annseptember        = 1041,
    annoctober          = 1042,
    annnovember         = 1043,
    anndecember         = 1044,
    sannjuly            = 1045,
    sannaugust          = 1046,
    sannseptember       = 1047,
    sannoctober         = 1048,
    sannnovember        = 1049,
    sanndecember        = 1050,
    unknown             =    0)
}

tifList <- function(){
  if(!exists(".tifList", envir = .tiEnv))
    tl <- setDefaultFrequencies(setup = TRUE)
  else{
    tl <- get(".tifList", envir = .tiEnv)
    if(is.na(tl["weekly"]))
      tl <- setDefaultFrequencies(setup = FALSE)
  }
  return(tl)
}

setDefaultFrequencies <- function(weekly     = "wmonday",
                                  biweekly   = "bw1wednesday",
                                  bimonthly  = "bimonthdecember",
                                  quarterly  = "qdecember",
                                  annual     = "anndecember",
                                  semiannual = "sanndecember",
                                  setup = FALSE){
  if(setup) tl <- initialTifList()
  else {
    if(!exists(".tifList", envir = .tiEnv))
      assign(".tifList", initialTifList(), envir = .tiEnv)
    tl <- get(".tifList", envir = .tiEnv)
  }
  if(!missing(weekly) || setup)
    tl["weekly"] <- tl[weekly]
  if(!missing(biweekly) || setup)
    tl["biweekly"] <- tl[biweekly]
  if(!missing(bimonthly) || setup)
    tl[c("bimonth", "bimonthly")] <- tl[bimonthly]
  if(!missing(quarterly) || setup)
    tl[c("quarterly", "q")] <- tl[quarterly]
  if(!missing(annual) || setup)
    tl[c("annual", "a")] <- tl[annual]
  if(!missing(semiannual) || setup)
    tl[c("sann", "semiannual")]  <- tl[semiannual]
  assign(".tifList", tl, envir = .tiEnv)
}


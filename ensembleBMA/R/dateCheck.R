`dateCheck` <-
function (YYYYMMDDHH)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
# require("chron")

 origin <- c(month = 1, day = 1, year = 2000)

 YYYYMMDDHH <- sapply(YYYYMMDDHH, as.character)
 chk <- rep(TRUE, length(YYYYMMDDHH))
 l <- sapply(YYYYMMDDHH, nchar)
 if (any(I <- (l > 10 | l == 9 | l < 8))) chk[I] <- FALSE
 dropHour <- l[!I] == 8
 
 year <- as.numeric(sapply( YYYYMMDDHH[!I], substring, first = 1, last = 4))
 month <- as.numeric(sapply( YYYYMMDDHH[!I], substring, first = 5, last = 6))
 day <- as.numeric(sapply( YYYYMMDDHH[!I], substring, first = 7, last = 8))
 julianDate0 <- julian( month, day, year, origin = origin)

 L <- length(YYYYMMDDHH[!I])
 hour <- rep( 0, L)
 hour[!dropHour] <- as.numeric(sapply((YYYYMMDDHH[!I])[!dropHour], 
                             substring, first = 9, last = 10))
 julianDate <- julianDate0 + hour/24
 ymdh <- rep( "YYYYMMDDHH", L)
 for (i in seq(along = ymdh)) {
    ymdh[i] <- julTOymdh( julianDate[i], origin=origin, dropHour=dropHour[i])
 }
 chk[!I] <- ymdh == YYYYMMDDHH[!I]
 chk
}


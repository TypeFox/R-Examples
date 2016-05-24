`ymdhTOjul` <-
function (YYYYMMDDHH, origin = c(month = 1, day = 1, year = 2000)) 
{
  #
  # copyright 2006-present, University of Washington. All rights reserved.
  # for terms of use, see the LICENSE file
  #             

# require("chron")

  YYYYMMDDHH <- sapply(YYYYMMDDHH, as.character)
  l <- unique(sapply(YYYYMMDDHH, nchar))
  if (length(l) > 1 || l > 10 || l == 9 || l < 8) {
    stop("input YYYYMMDDHH must be uniformly of length 8 or 10")
  }

  year <- as.numeric(sapply( YYYYMMDDHH, substring, first = 1, last = 4))
  month <- as.numeric(sapply( YYYYMMDDHH, substring, first = 5, last = 6))
  day <- as.numeric(sapply( YYYYMMDDHH, substring, first = 7, last = 8))
  julianDate0 <- julian( month, day, year, origin = origin)

  L <- length(YYYYMMDDHH)
  if (l == 8) {
    hour <- rep( 0, L)
  }
  else {
    hour <- as.numeric(sapply(YYYYMMDDHH, substring, first = 9, last = 10))
    if (any(hour >= 24)) stop("hour should be less than 24")
  }

  julianDate <- structure(julianDate0 + hour/24, origin = origin, nchar = l)

  ymdh <- julTOymdh( julianDate, dropHour = (l == 8), origin = origin)

  I <- ymdh == YYYYMMDDHH
  if (any(!I)) {
    print(YYYYMMDDHH[!I])
    stop("improper date(s)")
  }
  
 julianDate
}


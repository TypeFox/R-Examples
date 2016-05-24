nextBusinessDay <- function(x, holidays = NULL, goodFriday = F,
                            board = F, inaug = board){
  ## returns a business-day ti
  bDay <- ti(x, "business")
  weekend <- jul(x) != jul(bDay)
  z <- bDay - weekend + 1
  if(missing(holidays))
    holidays <- holidays(year(z), goodFriday = goodFriday,
                         board = board, inaug = inaug)
  hol <- match(ymd(z), holidays, nomatch = 0) > 0
  z[hol] <- z[hol] + 1
  z
}

previousBusinessDay <- function(x, holidays = NULL, goodFriday = F,
                                board = F, inaug = board){
  ## returns a business-day ti
  z <- ti(x, "business") - 1
  if(missing(holidays))
    holidays <- holidays(year(z), goodFriday = goodFriday,
                         board = board, inaug = inaug)
  hol <- match(ymd(z), holidays, nomatch = 0) > 0
  z[hol] <- z[hol] - 1
  z
}

isHoliday <- function(x, goodFriday = F, board = F, inaug = board, businessOnly = T){
  hols <- holidays(year(x), goodFriday, board, inaug, businessOnly)
  match(ymd(x), hols, nomatch = 0) > 0
}

isBusinessDay<- function(x, ...){
  dow <- dayOfWeek(x)
  return( dow > 1 & dow < 7 & !isHoliday(x, ...))
}

isGoodFriday <- function(x){
  hols <- goodFriday(year(x))
  match(ymd(x), hols, nomatch = 0) > 0
}

isEaster <- function(x){
  hols <- easter(year(x))
  match(ymd(x), hols, nomatch = 0) > 0
}

holidays <- function(years, goodFriday = F, board = F, inaug = board, businessOnly = T){
  ## Presidential Inauguration days are holidays in the DC area, but only if
  ## they fall on a weekday.  Inauguration day is usually Jan 20'th in the
  ## year following a year divisible by 4 (election years) except if Jan 20'th
  ## falls on a Sunday, in which case the Inauguration is held on the 21'st.
  ## Inauguration days will show if and only if inaug is TRUE.
  hols <- federalHolidays(years, board = board, businessOnly = businessOnly)
  if(goodFriday) hols <- sort(c(hols, goodFriday(years)))
  if(inaug){
    inaugDays <- inaugurationDay(years)
    if(length(inaugDays) > 0){
      inaugDays <- inaugDays[dayOfWeek(inaugDays) %in% 2:6]
      if(length(inaugDays) > 0)
        hols <- sort(c(hols, inaugDays))
    }
  }
  hols
}

federalHolidays <- function(years, board = F, businessOnly = T){
  ## returns yyyymmdd dates of federal holidays for given years.
  ##
  ## Federal law defines 10 holidays.  4 are set by date:
  ##
  ##   NewYears       January 1
  ##   Independence   July 4
  ##   Veterans       November 11
  ##   Christmas      December 25
  ##
  ## and the other 6 are set by day of the week and month:
  ##
  ##   MLK            third Monday of January
  ##   Presidents     third Monday of February
  ##   Memorial       last Monday of May
  ##   Labor          first Monday of September
  ##   Columbus       second Monday of October
  ##   Thanksgiving   fourth Thursday of November
  ##
  ## If one of the four fixed-date holidays falls on a Sunday, the federal
  ## holiday is celebrated the next day (Monday).  If it falls on a Saturday,
  ## the preceding day (Friday) is a holiday for the Federal Reserve Board,
  ## but not for the Reserve Banks and the banking system as a whole.
  ##
  ## If businessOnly is TRUE, drop the Saturday holidays. Note that this has
  ## no effect if board is TRUE, since that moves Saturday holidays to the
  ## preceding Friday

  z <- c(NewYearsDay(years),
         MLKingDay(years),
         GWBirthday(years),
         MemorialDay(years),
         IndependenceDay(years),
         LaborDay(years),
         ColumbusDay(years),
         VeteransDay(years),
         ThanksgivingDay(years),
         ChristmasDay(years))
  hols <- sort(z)
  weekday <- dayOfWeek(hols)
  if(any(weekday == 1))
    hols[weekday == 1] <- ymd(jul(hols[weekday == 1]) + 1)
  if(board && any(weekday == 7)){
    hols[weekday == 7] <- ymd(jul(hols[weekday == 7]) - 1)
    ## recompute weekday vector because Saturday holidays were moved to Friday
    weekday <- dayOfWeek(hols)
  }
  if(businessOnly){
    ## drop Saturdays (no need to drop Sundays since there aren't any)
    hols <- hols[weekday != 7]
  }
  hols
}

NewYearsDay <- function(years){
  years <- years[years > 1870]
  if(length(years) == 0) return(numeric(0))
  ans <- 10000*years + 101
  names(ans) <- rep("NewYears", length(ans))
  sort(ans)
}

MLKingDay <- function(years){
  years <- years[years > 1985]
  if(length(years) == 0) return(numeric(0))
  ans <- ymd(ti(10000*years + 115, "wmonday"))
  names(ans) <- rep("MLKing", length(ans))
  sort(ans)
}

GWBirthday <- function(years){
  pre1971 <- between(years, 1880, 1970)
  recent  <- years > 1970
  if(any(pre1971))
    ans <- 10000*years[pre1971] + 222
  else
    ans <- numeric(0)
  if(any(recent))
    ans <- c(ans, ymd(ti(10000*years[recent] + 215, "wmonday")))
  names(ans) <- rep("GWBirthday", length(ans))
  sort(ans)
}

MemorialDay <- function(years){
  pre1971 <- between(years, 1888, 1970)
  recent  <- years > 1970
  if(any(pre1971))
    ans <- 10000*years[pre1971] + 530
  else
    ans <- numeric(0)
  if(any(recent))
    ans <- c(ans, ymd(ti(10000*years[recent] + 601, "wmonday") - 1))
  names(ans) <- rep("Memorial", length(ans))
  sort(ans)
}

IndependenceDay <- function(years){
  years <- years[years >= 1870]
  if(length(years) == 0) return(numeric(0))
  ans <- 10000*years + 704
  names(ans) <- rep("Independence", length(ans))
  sort(ans)
}

LaborDay <- function(years){
  years <- years[years >= 1894]
  if(length(years) == 0) return(numeric(0))
  ans <- ymd(ti(years*10000 +  901, "wmonday"))
  names(ans) <- rep("Labor", length(ans))
  sort(ans)
}

ColumbusDay <- function(years){
  pre1971 <- between(years, 1934, 1970)
  recent  <- years > 1970
  if(any(pre1971))
    ans <- 10000*years[pre1971] + 1012
  else
    ans <- numeric(0)
  if(any(recent))
    ans <- c(ans, ymd(ti(10000*years[recent] + 1008, "wmonday")))
  names(ans) <- rep("Columbus", length(ans))
  sort(ans)
}

VeteransDay <- function(years){
  mondayYears <- between(years, 1971, 1977)
  nov11Years <- between(years, 1938, 1970) | years > 1977
  if(any(mondayYears))
    ans <- ymd(ti(10000*years[mondayYears] + 1022, "wmonday"))
  else
    ans <- numeric(0)
  if(any(nov11Years))
    ans <- c(ans, 10000*years[nov11Years] + 1111)
  names(ans) <- rep("Veterans", length(ans))
  sort(ans)
}

ThanksgivingDay <- function(years){
  pre1939Years <- years[between(years, 1863, 1938)]
  year1939 <- years[years == 1939]
  year1940 <- years[years == 1940]
  recentYears <- years[years > 1940]
  if(length(pre1939Years) > 0){
    ## last Thursday of November, can't use ti function since years
    ## may be earlier than 1900 
    novEnds <- 10000*pre1939Years + 1130
    ans <- novEnds - ((dayOfWeek(novEnds) + 2) %% 7)
  }
  else
    ans <- numeric(0)
  if(length(year1939) > 0)
    ans <- c(ans, rep(19391123, length(year1939)))
  if(length(year1940) > 0)
    ans <- c(ans, rep(19401121, length(year1940)))
  if(length(recentYears) > 0)
    ans <- c(ans, ymd(ti(10000*recentYears + 1122, "wthursday")))
  names(ans) <- rep("Thanksgiving", length(ans))
  sort(ans)
}

ChristmasDay <- function(years){
  years <- years[years >= 1870]
  if(length(years) == 0) return(numeric(0))
  ans <- 10000*years + 1225
  names(ans) <- rep("Christmas", length(ans))
  sort(ans)
}

goodFriday <- function(years){
  ## yyyymmdd dates of Good Friday for given years
  z <- ymd(jul(easter(years)) - 2)
  names(z) <- rep("GoodFriday", length(z))
  z
}

inaugurationDay <- function(years){
  ## yyyymmdd dates of Inauguration Day for given years
  inaugDates <- c(17890430,
                  10000*(1793 + 4*(0:35))  + 304,
                  10000*(1937 + 4*(0:115)) + 120)
  jInaug <- jul(as.Date(as.character(inaugDates), format = "%Y%m%d"))
  sunday <- unclass(jInaug) %% 7 == 6
  jInaug[sunday] <- jInaug[sunday] + 1
  keep <- which(year(jInaug) %in% years)
  if(length(which) == 0)
    return(numeric(0))
  else {
    z <- sort(ymd(jInaug[keep]))
    names(z) <- rep("Inauguration", length(z))
    return(z)
  }
}

easter <- function(years){
  ## yyyymmdd dates of Easter for supplied 4 digit years
  G <- years %% 19
  C <- years %/% 100
  H <- (C - (C %/% 4) - ((8*C + 13) %/% 25) + 19*G + 15) %% 30
  I <- H - (H %/% 28) * (1 - (H %/% 28)*(29 %/% (H + 1))*((21 - G) %/% 11))
  J <- (years + (years %/% 4) + I + 2 - C + (C %/% 4)) %% 7
  L <-  I - J
  month <- 3 + (L + 40) %/% 44
  day <- L + 28 - 31*(month %/% 4)
  10000*years + 100*month + day
}

holidaysBetween <- function(startTi, endTi, goodFriday = F,
                            board = F,  inaug = board, businessOnly = T){
  startTi <- tiDaily(startTi)
  endTi   <- tiDaily(endTi)
  years <- year(startTi):year(endTi)
  hols <- holidays(years, goodFriday, board, inaug, businessOnly)
  hols[between(hols, ymd(startTi), ymd(endTi))]
}

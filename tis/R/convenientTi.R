## The argument xTi in the functions below is a ti (time index).

## dates used in various reports, releases, etc...
today <- function(tif = "daily"){
  ## today's date as a daily ti
  ti(Sys.Date(), tif = tif)
}

dayOfPeriod <- function(xTi = today(), tif = NULL){
  if(is.null(tif)) stop("Bad or missing tif")
  julianDates <- floor(jul(xTi))
  prevPeriodEnds <- jul(ti(julianDates, tif) - 1)
  julianDates - prevPeriodEnds
}

dayOfWeek <- function(xTi = today()){
  ## 1 = Sunday, 7 = Saturday
  (floor(unclass(jul(xTi) + 1)) %% 7) + 1
}

dayOfMonth <- function(xTi = today()){
  floor(ymd(xTi)) %% 100
}

dayOfYear <- function(xTi = today()){
  floor(jul(xTi)) + 1 - jul(10000*year(xTi) + 101)
}

tiDaily <- function(xTi, offset = 1){
  ## returns daily ti(s) 
  ti(jul(xTi, offset), "daily")
}

tiBusiness <- function(xTi, offset = 1){
  ## returns business ti(s) 
  ti(jul(xTi, offset), "business")
}

firstDayOf <- function(xTi)
  ti(jul(xTi - 1) + 1, "daily")

lastDayOf <- function(xTi)
  ti(jul(xTi), "daily")

firstBusinessDayOf <- function(xTi)
  ti(jul(xTi - 1) + 1, "business")

lastBusinessDayOf <- function(xTi)
  firstBusinessDayOf(xTi + 1) - 1

lastBusinessDayOfMonth <- function(xTi){
  ## return ti for last business day of the month of xTi
  lastBusinessDayOf(ti(xTi, "monthly"))
}

firstBusinessDayOfMonth <- function(xTi){
  ## return ti for first business day of the month of xTi
  firstBusinessDayOf(ti(xTi, "monthly"))
}

currentMonthDay <- function(xTi, daynum){
  ## the next upcoming daynum'th of the month
  tiDaily(currentMonth(tiDaily(xTi) - daynum)) + daynum  
}

latestMonthDay <- function(xTi, daynum){
  ## the most recent daynum'th of the month
  tiDaily(latestMonth(tiDaily(xTi) - daynum)) + daynum  
}

## currentwhatever(xTi) 
currentMonday    <- function(xTi = today()) tiDaily(ti(xTi, "wmonday"))
currentTuesday   <- function(xTi = today()) tiDaily(ti(xTi, "wtuesday"))
currentWednesday <- function(xTi = today()) tiDaily(ti(xTi, "wwednesday"))
currentThursday  <- function(xTi = today()) tiDaily(ti(xTi, "wthursday"))
currentFriday    <- function(xTi = today()) tiDaily(ti(xTi, "wfriday"))
currentSaturday  <- function(xTi = today()) tiDaily(ti(xTi, "wsaturday"))
currentSunday    <- function(xTi = today()) tiDaily(ti(xTi, "wsunday"))
currentWeek      <- function(xTi = today()) ti(xTi, "weekly")
currentMonth     <- function(xTi = today()) ti(xTi, "monthly")   
currentQuarter   <- function(xTi = today()) ti(xTi, "quarterly") 
currentHalf      <- function(xTi = today()) ti(xTi, "semiannual") 
currentYear      <- function(xTi = today()) ti(xTi, "annual")    
currentQ4        <- function(xTi = today()) ti(jul(currentYear(xTi)), "quarterly")
currentQMonth    <- function(xTi = today()) currentMonth(currentQuarter(xTi))

currentJanuary   <- function(xTi = today()) ti(ti(xTi, "annjanuary"),   "monthly")
currentFebruary  <- function(xTi = today()) ti(ti(xTi, "annfebruary"),  "monthly")
currentMarch     <- function(xTi = today()) ti(ti(xTi, "annmarch"),     "monthly")
currentApril     <- function(xTi = today()) ti(ti(xTi, "annapril"),     "monthly")
currentMay       <- function(xTi = today()) ti(ti(xTi, "annmay"),       "monthly")
currentJune      <- function(xTi = today()) ti(ti(xTi, "annjune"),      "monthly")
currentJuly      <- function(xTi = today()) ti(ti(xTi, "annjuly"),      "monthly")
currentAugust    <- function(xTi = today()) ti(ti(xTi, "annaugust"),    "monthly")
currentSeptember <- function(xTi = today()) ti(ti(xTi, "annseptember"), "monthly")
currentOctober   <- function(xTi = today()) ti(ti(xTi, "annoctober"),   "monthly")
currentNovember  <- function(xTi = today()) ti(ti(xTi, "annnovember"),  "monthly")
currentDecember  <- function(xTi = today()) ti(ti(xTi, "anndecember"),  "monthly")

## latestwhatever(xTi)
latestMonday    <- function(xTi = today()) tiDaily(ti(jul(xTi) - 6, "wmonday"))
latestTuesday   <- function(xTi = today()) tiDaily(ti(jul(xTi) - 6, "wtuesday"))
latestWednesday <- function(xTi = today()) tiDaily(ti(jul(xTi) - 6, "wwednesday"))
latestThursday  <- function(xTi = today()) tiDaily(ti(jul(xTi) - 6, "wthursday"))
latestFriday    <- function(xTi = today()) tiDaily(ti(jul(xTi) - 6, "wfriday"))
latestSaturday  <- function(xTi = today()) tiDaily(ti(jul(xTi) - 6, "wsaturday"))
latestSunday    <- function(xTi = today()) tiDaily(ti(jul(xTi) - 6, "wsunday"))
latestWeek      <- function(xTi = today()) ti(jul(xTi) + 1, "weekly") - 1
latestMonth     <- function(xTi = today()) ti(jul(xTi) + 1, "monthly") - 1
latestQuarter   <- function(xTi = today()) ti(jul(xTi) + 1, "quarterly") - 1
latestHalf      <- function(xTi = today()) ti(jul(xTi) + 1, "semiannual") - 1
latestYear      <- function(xTi = today()) ti(jul(xTi) + 1, "annual") - 1
latestQ4        <- function(xTi = today()) ti(jul(latestYear(xTi)), "quarterly")

latestJanuary    <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annjanuary")   - 1, "monthly")
latestFebruary   <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annfebruary")  - 1, "monthly")
latestMarch      <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annmarch")     - 1, "monthly")
latestApril      <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annapril")     - 1, "monthly")
latestMay        <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annmay")       - 1, "monthly")
latestJune       <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annjune")      - 1, "monthly")
latestJuly       <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annjuly")      - 1, "monthly")
latestAugust     <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annaugust")    - 1, "monthly")
latestSeptember  <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annseptember") - 1, "monthly")
latestOctober    <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annoctober")   - 1, "monthly")
latestNovember   <- function(xTi = today()) ti(ti(jul(xTi) + 1, "annnovember")  - 1, "monthly")
latestDecember   <- function(xTi = today()) ti(ti(jul(xTi) + 1, "anndecember")  - 1, "monthly")


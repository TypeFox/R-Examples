library(PCICt)
library(RUnit)

test.values <- c("1850-01-01", "2012-04-01 00:10:00", "2012-02-28", "1599-01-01", "1582-10-04", "1582-10-03", "900-02-28", "900-03-01", "200-03-01", "1-01-01", "0-01-01")

## Tests as.pcict functions
PCICt.test.as.PCICt <- function() {
  test.calendars <- c("360", "365", "gregorian", "proleptic_gregorian")
  cal.360.only <- c("2011-02-30")
  gregorian.not.365 <- c("2012-02-29")
  gregorian.not.360 <- c("2011-01-31")

  ## Check that valid input produces valid output for all calendar types
  for(cal in test.calendars) {
    for(i in seq_along(test.values)) {
      checkEquals(as.character(as.PCICt(test.values[i], cal=cal)), test.values[i])
      
      ## Check that as.PCICt.POSIXlt works as expected
      checkEquals(as.character(as.PCICt(as.POSIXlt(test.values[i]), cal=cal)), test.values[i])
      
      ## Check that as.PCICt.POSIXct works as expected
      ##checkEquals(as.character(as.PCICt(as.POSIXct(test.values[i]), cal=cal)), test.values[i])
    }
  }

  ## Check factor case
  checkEquals(as.PCICt(as.factor(gregorian.not.365), cal="360"), as.PCICt(gregorian.not.365, cal="360"))
  
  ## Check that invalid input triggers errors
  checkException(as.PCICt(cal.360.only, cal="365"))
  checkException(as.PCICt(cal.360.only, cal="gregorian"))

  ## Should fail but doesn't; not sure how to hack it other than to pass in month lengths to the so-called "360" funcs.
  ##checkException(as.PCICt(gregorian.not.365, cal="365"))

  checkException(as.PCICt(gregorian.not.360, cal="360"))
  checkException(as.PCICt("your mom", cal="360"))
  checkException(as.PCICt("your mom", cal="365"))

  ## Check that NAs are passed through properly
  checkTrue(is.na(as.PCICt(NA, cal="360")))
  checkTrue(is.na(as.PCICt(NA, cal="365")))

  ## Check numeric input...
  checkEquals(as.PCICt(86400, origin="1961-01-01", cal="360"), as.PCICt("1961-01-02", cal="360"))
}

## Tests subset operators
PCICt.test.subset <- function() {
  dat <- as.PCICt(test.values, cal="365")
  dat2 <- dat[2:4]
  dat3 <- dat
  dat3[1:2] <- dat[2:3]
  
  ## Check that subset operator works as expected
  checkEquals(as.character(dat[2:4]), as.character(dat)[2:4])

  ## Check that subset operator preserves attributes properly
  checkEquals(attr(dat2, "cal"), attr(dat, "cal"))
  checkEquals(attr(dat3, "cal"), attr(dat, "cal"))

  ## Check assignment
  checkEquals(dat3[1:4], dat[c(2, 3, 3, 4)])

  ## Check exception handling
  dat.360 <- as.PCICt(test.values, cal="360")
  checkException(dat3[1:3] <- dat.360[2:4])
}

## Tests arithmetic operators (+, -)
PCICt.test.operators <- function() {
  test.365 <- as.PCICt(test.values, cal="365")
  test.greg <- as.PCICt(test.values, cal="gregorian")
  
  ## Check edge cases for each calendar for addition (isn't February a great month?)
  checkEquals(as.character(as.PCICt("2012-02-28", cal="365") + 86400), "2012-03-01")
  checkEquals(as.character(as.PCICt("2011-02-28", cal="360") + 2 * 86400), "2011-02-30")
  checkEquals(as.character(as.PCICt("2012-02-28", cal="gregorian") + 86400), "2012-02-29")

  ## FIXME: CHECK THAT SUBTRACTION GIVES A DIFFTIME OBJECT
  checkException(test.365 - test.greg)
  checkException(test.365 - test.greg)
  
  ## Check edge cases for each calendar for subtraction
  checkEquals(as.numeric(as.PCICt("2012-03-01", cal="365") - as.PCICt("2012-02-28", cal="365"), units="secs"), as.numeric(as.difftime("24:00:00"), units="secs"))
  checkEquals(as.numeric(as.PCICt("2011-03-01", cal="360") - as.PCICt("2011-02-30", cal="360"), units="secs"), as.numeric(as.difftime("24:00:00"), units="secs"))
  checkEquals(as.numeric(as.PCICt("2012-03-01", cal="gregorian") - as.PCICt("2012-02-29", cal="gregorian"), units="secs"), as.numeric(as.difftime("24:00:00"), units="secs"))
}

## Tests the truncate function
PCICt.test.trunc <- function() {
  date.360 <- as.PCICt("2012-02-28 12:34:56.7", cal="360")
  date.365 <- as.PCICt("2012-02-28 12:34:56.7", cal="365")
  date.gregorian <- as.PCICt("2012-02-28 12:34:56.7", cal="gregorian")

  results <- c("2012-02-28 12:34:56", "2012-02-28 12:34:00", "2012-02-28 12:00:00", "2012-02-28")
  round.to <- c("secs", "mins", "hours", "days")
  
  for(i in 1:length(results)) {
    checkEquals(as.character(trunc(date.360, round.to[i])), results[i])
    checkEquals(as.character(trunc(date.365, round.to[i])), results[i])
    checkEquals(as.character(trunc(date.gregorian, round.to[i])), results[i])
  }
}

## Tests the round function
PCICt.test.round <- function() {
  date.360 <- as.PCICt("2012-02-28 12:34:56.7", cal="360")
  date.365 <- as.PCICt("2012-02-28 12:34:56.7", cal="365")
  date.gregorian <- as.PCICt("2012-02-28 12:34:56.7", cal="gregorian")

  results <- c("2012-02-28 12:34:57", "2012-02-28 12:35:00", "2012-02-28 13:00:00")
  results.days <- c("2012-02-29", "2012-03-01", "2012-02-29")
  round.to <- c("secs", "mins", "hours")
  
  for(i in 1:length(results)) {
    checkEquals(as.character(round(date.360, round.to[i])), results[i])
    checkEquals(as.character(round(date.365, round.to[i])), results[i])
    checkEquals(as.character(round(date.gregorian, round.to[i])), results[i])
  }
  checkEquals(as.character(round(date.360, "days")), results.days[1])
  checkEquals(as.character(round(date.365, "days")), results.days[2])
  checkEquals(as.character(round(date.gregorian, "days")), results.days[3])
}

## Tests the cut function
PCICt.test.cut <- function() {
}

PCICt.test.summary <- function() {
}

PCICt.test.format <- function() {
}

PCICt.test.as.POSIXlt <- function() {
}

PCICt.test.as.POSIXct <- function() {
}

PCICt.test.julian <- function() {
}

PCICt.test.range <- function() {
}

PCICt.test.pretty <- function() {
}

PCICt.test.c <- function() {
}

PCICt.test.seq <- function() {
  ## Check that sequences work as expected
  checkEquals(seq(as.PCICt("2011-02-30", cal="360"), as.PCICt("2012-04-01 00:10:00", cal="360"), by="years"), as.PCICt(c("2011-02-30", "2012-02-30"), cal="360"))
  checkEquals(seq(as.PCICt("2011-02-28", cal="365"), as.PCICt("2012-04-01 00:10:00", cal="365"), by="years"), as.PCICt(c("2011-02-28", "2012-02-28"), cal="365"))
  checkEquals(seq(as.PCICt("2011-02-28", cal="gregorian"), as.PCICt("2012-04-01 00:10:00", cal="gregorian"), by="years"), as.PCICt(c("2011-02-28", "2012-02-28"), cal="gregorian"))

  checkEquals(seq(as.PCICt("2011-02-30", cal="360"), as.PCICt("2011-03-03 00:10:00", cal="360"), by="days"), as.PCICt(c("2011-02-30", "2011-03-01", "2011-03-02", "2011-03-03"), cal="360"))
  checkEquals(seq(as.PCICt("2012-02-28", cal="365"), as.PCICt("2012-03-03 00:10:00", cal="365"), by="days"), as.PCICt(c("2012-02-28", "2012-03-01", "2012-03-02", "2012-03-03"), cal="365"))
  checkEquals(seq(as.PCICt("2012-02-28", cal="gregorian"), as.PCICt("2012-03-03 00:10:00", cal="gregorian"), by="days"), as.PCICt(c("2012-02-28", "2012-02-29", "2012-03-01", "2012-03-02", "2012-03-03"), cal="gregorian"))

}

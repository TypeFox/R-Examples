library(climdex.pcic)
library(RUnit)

climdex.pcic.test.quantile <- function() {
  ## Check a basic but fundamental case.
  dat <- c(5, 8, 1, 2, 10, 3, 9, 4, 7, 6)
  qtiles <- c(0, 0.1, 0.5, 0.9, 1.0)
  
  checkEquals(as.numeric(quantile(dat, qtiles, type=8)), climdex.quantile(dat, qtiles))

  ## Test the n=1 case
  checkEquals(as.numeric(quantile(1, qtiles, type=8)), climdex.quantile(1, qtiles))

  ## Test the n=0 case
  checkEquals(as.numeric(quantile(c(), qtiles, type=8)), climdex.quantile(c(), qtiles))
}

climdex.pcic.test.select.blocks <- function() {
  ## Check a case that should prove problematic.
  dat1 <- c(rep(TRUE, 4), FALSE, rep(TRUE, 3), rep(FALSE, 2), rep(TRUE, 5))
  result1 <- c(rep(TRUE, 4), rep(FALSE, 4), rep(FALSE, 2), rep(TRUE, 5))
  dat2 <- c(FALSE, rep(TRUE, 4), FALSE, rep(TRUE, 3), rep(FALSE, 2), rep(TRUE, 5), FALSE)
  result2 <- c(FALSE, rep(TRUE, 4), rep(FALSE, 4), rep(FALSE, 2), rep(TRUE, 5), FALSE)
  checkEquals(result1, select.blocks.gt.length(dat1, 3))
  checkEquals(result2, select.blocks.gt.length(dat2, 3))
  checkEquals(dat1, select.blocks.gt.length(dat1, 0))
}

climdex.pcic.test.total.precip.op.threshold <- function() {
  prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
  phony.date.factor <- factor(rep(1:2, each=5))
  checkEquals(as.numeric(total.precip.op.threshold(prec.dat, phony.date.factor, 1, ">=")), c(10.5, 11))
}

climdex.pcic.test.spell.length.max <- function() {
  prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
  phony.date.factor <- factor(rep(1:2, each=5))
  cwd.result <- c(0, 5)
  altcwd.result <- c(4, 2)
  
  ## With spells spanning years...
  cwd <- as.numeric(spell.length.max(prec.dat, phony.date.factor, 1, ">=", TRUE))
  
  ## Without spells spanning years...
  altcwd <- as.numeric(spell.length.max(prec.dat, phony.date.factor, 1, ">=", FALSE))

  checkEquals(cwd.result, cwd)
  checkEquals(altcwd.result, altcwd)
}

climdex.pcic.test.simple.precipitation.intensity.index <- function() {
  prec.dat <- c(0.1, 3.0, 4.3, 0.9, 1.3, 6.0, 0, 0, 4.0, 1)
  phony.date.factor <- factor(rep(1:2, each=5))
  sdii <- as.numeric(simple.precipitation.intensity.index(prec.dat, phony.date.factor))
  sdii.result <- c(43/15, 11/3)
  checkEquals(sdii.result, sdii)
}

climdex.pcic.test.nday.consec.prec.max <- function() {
  ##              1    1    1    1    1    1  1  1  1  1   1  1   2  2  2  2  2  2   2  2  2  2  2   2
  prec.dat <- c(0.1, 3.0, 4.3, 0.9, 0.7, 1.7, 8, 1, 0, 0, 30, 0, 10, 8, 0, 0, 0, 1, 20, 0, 0, 0, 0, 50)
  phony.date.factor <- factor(rep(1:2, each=12))

  result.1day <- c(30, 50)
  result.5day.center <- c(48, 50)
  result.5day.end <- c(39, 50)
  result.6day <- c(48, 70)

  checkEquals(result.1day, as.numeric(nday.consec.prec.max(prec.dat, phony.date.factor, 1)))
  checkEquals(result.5day.center, as.numeric(nday.consec.prec.max(prec.dat, phony.date.factor, 5)))
  checkEquals(result.5day.end, as.numeric(nday.consec.prec.max(prec.dat, phony.date.factor, 5, TRUE)))
  checkEquals(result.6day, as.numeric(nday.consec.prec.max(prec.dat, phony.date.factor, 6)))
}

climdex.pcic.test.threshold.exceedance.duration.index <- function() {
  prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
  phony.date.factor <- factor(rep(1:2, each=5))
  phony.jdays <- rep(1:5, 2)
  thresholds <- c(0.5, 1, 1.5, 1.5, 0.5)
  op <- ">="
  min.length <- 2
  
  ## With spells spanning years...
  alttedi <- threshold.exceedance.duration.index(prec.dat, phony.date.factor, phony.jdays, thresholds, op, min.length, TRUE, 1)

  ## Without spells spanning years...
  tedi <- threshold.exceedance.duration.index(prec.dat, phony.date.factor, phony.jdays, thresholds, op, min.length, FALSE, 1)

  checkEquals(c(4, 3), as.numeric(alttedi))
  checkEquals(c(4, 2), as.numeric(tedi))
}

climdex.pcic.test.get.series.lengths.at.ends <- function() {
  ## Get lengths of sequences of TRUE values in a sequence
  series.lengths <- get.series.lengths.at.ends(c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))

  checkEquals(c(0, 0, 3, 0, 1, 0, 0, 0, 0, 4, 0), series.lengths)
}

climdex.pcic.test.number.days.op.threshold <- function() {
  data <- c(1, 3, 4, 0, 1, 2, 3, 4)
  date.factor <- factor(rep(1:2, each=4))
  op <- ">="
  threshold <- 2
  checkEquals(c(2, 3), as.numeric(number.days.op.threshold(data, date.factor, threshold, op)))
}

climdex.pcic.test.percent.days.op.threshold <- function() {
  is.x86 <- !(climdex.quantile(c(0, 1, 2), 0.3) == quantile(c(0, 1, 2), 0.3, type=8))
  
  ## Parse the dates into PCICt.
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  
  ## Load the data in.
  ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP, ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION, tmax.dates, tmin.dates, prec.dates, base.range=c(1981, 1990))
 
  ## Compute monthly tx90p on example data.
  tx90p <- percent.days.op.threshold(ci@data$tmax, ci@dates, ci@jdays, ci@date.factors$monthly, ci@quantiles$tmax$outbase$q90, ci@quantiles$tmax$inbase$q90, ci@base.range, ">", ci@max.missing.days['monthly']) * ci@namasks$monthly$tmax

  ## Check this against the known-valid data.
  if(!is.x86)
    checkEquals(ec.1018935.tx90p, tx90p)
}

climdex.pcic.test.growing.season.length <- function() {
  ## Parse the dates into PCICt.
 tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
 tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
 prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year", "jday")]), format="%Y %j", cal="gregorian")
 
 ## Load the data in.
 ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP, ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION, tmax.dates, tmin.dates, prec.dates, base.range=c(1981, 1990))
 
 ## Create an annual timeseries of the growing season length in days.
 gsl <- growing.season.length(ci@data$tavg, ci@date.factors$annual, ci@dates, ci@northern.hemisphere, gsl.mode="GSL") * ci@namasks$annual$tavg

 ## Check this against known-valid data.
 checkEquals(ec.1018935.gsl, gsl)

 ## Now for specific checks of the various modes...
 test.threshold <- 3
 lv <- test.threshold - 1
 gv <- test.threshold + 1
 ls <- 365
 num.days.set <- c(1, 3, 5, 7, 9)
 names(num.days.set) <- num.days.set
 modes <- c("GSL", "GSL_first", "GSL_max", "GSL_sum")
 test.dates <- seq(as.PCICt("1961-01-01", cal=ls), as.PCICt("1961-12-31", cal=ls), by="days")
 test.dates.factor <- factor(format(test.dates, "%Y"))
 td <- list(t1=rep(lv, ls), t2=rep(gv, ls))
 td$t1[c(1:3, 8:15, 20:32, 40:70, 80:240, 250:270, 280:285)] <- gv
 td$t2[c(7:12, 40:80, 352:358)] <- lv
 
 res <- sapply(num.days.set, function(num.days) { sapply(td, function(dat) { sapply(modes, function(mode) { growing.season.length(dat, test.dates.factor, test.dates, TRUE, num.days, test.threshold, mode) }, simplify="array") }, simplify="array") }, simplify="array")
 expected.res <- structure(c(240, 3, 161, 243, 351, 6, 271, 311,
                             240, 3, 161, 243, 351, 6, 271, 311,
                             233, 25, 161, 244, 351, 6, 271, 311,
                             233, 25, 161, 238, 339, 27, 271, 305,
                             221, 51, 161, 233, 353, 27, 271, 298),
                           .Dim = c(4L, 2L, 5L),
                           .Dimnames = list(c("GSL.1961", "GSL_first.1961", "GSL_max.1961", "GSL_sum.1961"), c("t1", "t2"), c("1", "3", "5", "7", "9")))
}

climdex.pcic.test.bootstrap.quantiles <- function() {
  ## Calculates quantiles for days with an n-day window, with or without year replacement (leave out rm.yr.idx and dup.yr.idx to not do year replacement).
  get.quantile.for.day <- function(ci, win.size, q, var, day, rm.yr.idx, dup.yr.idx, dpy=365) {
    dat <- ci@data[[var]][climdex.pcic:::get.bootstrap.set(ci@dates, ci@base.range, win.size)]
    dim(dat) <- c(day=dpy, yr=length(dat) / dpy)
    half.win <- floor(win.size / 2);
    day.idx <- ((-half.win):half.win + (day - 1) + dpy) %% dpy + 1
    num.yrs <- dim(dat)["yr"]
    yr.idx <- if(missing(rm.yr.idx) || missing(dup.yr.idx)) 1:num.yrs else c((1:num.yrs)[-rm.yr.idx], dup.yr.idx)
    return(quantile(dat[day.idx, yr.idx], probs=q, type=8, na.rm=TRUE))
  }
  
  ## Checks that climdex.pcic in-base quantiles exactly match the (mostly) independent implementation above.
  check.inbase.quantiles.for.day <- function(ci, pctile, var, day) {
    climdex.pcic.data <- ci@quantiles[[var]]$inbase[[paste("q", pctile, sep="")]][day,,]
    nyr <- dim(climdex.pcic.data)[1]
    comparison.data <- t(sapply(1:nyr, function(x) { sapply((1:nyr)[-x], function(y) { get.quantile.for.day(ci, 5, pctile / 100, var, day, x, y) }) }))
    dimnames(comparison.data) <- NULL
    names(comparison.data) <- NULL
    checkEquals(climdex.pcic.data, comparison.data)
  }
  
  ## Checks that climdex.pcic out-of-base quantiles exactly match the (mostly) independent implementation above.
  check.outbase.quantiles <- function(ci, pctile, var) {
    climdex.pcic.data <- ci@quantiles[[var]]$outbase[[paste("q", pctile, sep="")]]
    comparison.data <- sapply(1:length(climdex.pcic.data), function(x) { get.quantile.for.day(ci, 5, pctile / 100, var, x) })
    dimnames(comparison.data) <- NULL
    names(comparison.data) <- NULL
    checkEquals(climdex.pcic.data, comparison.data)
  }
  
  ## Read in data, create climdexInput data structure.
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  
  ## Load the data in.
  ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP, ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION, tmax.dates, tmin.dates, prec.dates, base.range=c(1981, 1990))
  
  ## Check quantiles for 10/90 for tmin/tmax
  for(var in c("tmin", "tmax")) {
    for(pctile in c(10, 90)) {
      for(day in c(1, 2, 185, 364, 365)) {
        cat(paste("Checking in base for ", var, ", ", pctile, "th percentile, day ", day, "\n", sep=""))
        check.inbase.quantiles.for.day(ci, pctile, var, day)
      }
      cat(paste("Checking out base for ", var, ", ", pctile, "th percentile\n", sep=""))
      check.outbase.quantiles(ci, pctile, var)
    }
  }
}

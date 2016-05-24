library(climdex.pcic)
library(RUnit)

climdex.pcic.test.intake.routines <- function() {
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmax.dat <- ec.1018935.tmax$MAX_TEMP
  tmin.dat <- ec.1018935.tmin$MIN_TEMP
  prec.dat <- ec.1018935.prec$ONE_DAY_PRECIPITATION
  all.indices <- c('fd', 'su', 'id', 'tr', 'gsl', 'txx', 'tnx', 'txn', 'tnn', 'tn10p', 'tx10p', 'tn90p', 'tx90p', 'wsdi', 'csdi',
                   'dtr', 'rx1day', 'rx5day', 'sdii', 'r10mm', 'r20mm', 'rnnmm', 'cdd', 'cwd', 'r95ptot', 'r99ptot', 'prcptot')

  threshold.indices.to.exclude.on.x86 <- c("tx10p", "tx90p", "wsdi", "tn10p", "tn90p", "csdi", "r95ptot", "r99ptot")
  is.x86 <- !(climdex.quantile(c(0, 1, 2), 0.3) == quantile(c(0, 1, 2), 0.3, type=8))
  
  for(i in 1:7) {
    include.tmax <- i %% 2
    include.tmin <- floor(i / 2) %% 2
    include.prec <- floor(i / 4) %% 2
    ci <- climdexInput.raw(if(include.tmax) tmax.dat else NULL,
                           if(include.tmin) tmin.dat else NULL,
                           if(include.prec) prec.dat else NULL,
                           if(include.tmax) tmax.dates else NULL,
                           if(include.tmin) tmin.dates else NULL,
                           if(include.prec) prec.dates else NULL,
                           base.range=c(1981, 1990))

    outbase.thresholds <- get.outofbase.quantiles(if(include.tmax) tmax.dat else NULL,
                           if(include.tmin) tmin.dat else NULL,
                           if(include.prec) prec.dat else NULL,
                           if(include.tmax) tmax.dates else NULL,
                           if(include.tmin) tmin.dates else NULL,
                           if(include.prec) prec.dates else NULL,
                           base.range=c(1981, 1990))

    ci.csv <- climdexInput.csv(if(include.tmax) "1018935_MAX_TEMP.csv" else NULL,
                               if(include.tmin) "1018935_MIN_TEMP.csv" else NULL,
                               if(include.prec) "1018935_ONE_DAY_PRECIPITATION.csv" else NULL,
                               data.columns=list(tmax="MAX_TEMP", tmin="MIN_TEMP", prec="ONE_DAY_PRECIPITATION"),
                               base.range=c(1981, 1990))

    indices.to.check.equals <- climdex.get.available.indices(ci, function.names=FALSE)
    indices.to.check.error <- all.indices[!(all.indices %in% indices.to.check.equals)]

    ## 80-bit vs 64-bit SSE floating point differences cause comparisons to fail, thanks
    ## to threshold-based comparisons amplifying the effect.
    if(is.x86)
      indices.to.check.equals <- indices.to.check.equals[!(indices.to.check.equals %in% threshold.indices.to.exclude.on.x86)]
    
    print("Checking success...")
    for(index in indices.to.check.equals) {
      print(index)
      fun <- match.fun(paste('climdex', index, sep="."))
      valid.result <- get(paste('ec.1018935', index, sep="."))
      checkEquals(valid.result, fun(ci), paste(index, "didn't match expected result: include.tmax is", include.tmax, ", include.tmin is ", include.tmin, ", include.prec is ", include.prec))
      checkEquals(valid.result, fun(ci.csv), paste(index, "didn't match expected result with CSV input: include.tmax is", include.tmax, ", include.tmin is ", include.tmin, ", include.prec is ", include.prec))
    }

    print("Checking error...")
    for(index in indices.to.check.error) {
      print(index)
      fun <- match.fun(paste('climdex', index, sep="."))
      checkException(fun(ci), paste(index, "didn't produce error when it should have: include.tmax is", include.tmax, ", include.tmin is ", include.tmin, ", include.prec is ", include.prec))
      checkException(fun(ci.csv), paste(index, "didn't produce error when it should have with CSV input: include.tmax is", include.tmax, ", include.tmin is ", include.tmin, ", include.prec is ", include.prec))
    }
  }
}

climdex.pcic.test.thresholds.edge.cases <- function() {
  dat1 <- c(rep(1:10, 10), rep(NA, 53), rep(1:4, 53))
  dat1.dates <- seq(from=as.PCICt("1970-01-01", cal=365), to=as.PCICt("1970-12-31", cal=365), by="days")
  dat2 <- rep(c(rep(1:10, 10), rep(NA, 53), rep(1:4, 53)), 2)
  dat2.dates <- seq(from=as.PCICt("1970-01-01", cal=365), to=as.PCICt("1971-12-31", cal=365), by="days")
  dat3 <- rep(c(rep(1:38, 4), rep(NA, 5), rep(1:4, 52)), 3)
  dat3.dates <- seq(from=as.PCICt("1970-01-01", cal=365), to=as.PCICt("1972-12-31", cal=365), by="days")
  dat3[c(155 + 365 * 2 + c(0:1))] <- 12
  
  ci <- climdexInput.raw(tmax=dat1, tmin=dat1, tmax.dates=dat1.dates, tmin.dates=dat1.dates, base.range=c(1970, 1970))
  ci2 <- climdexInput.raw(tmax=dat2, tmin=dat2, tmax.dates=dat2.dates, tmin.dates=dat2.dates, base.range=c(1970, 1971))
  ci3 <- climdexInput.raw(tmax=dat3, tmin=dat3, tmax.dates=dat3.dates, tmin.dates=dat3.dates, base.range=c(1970, 1971))

  checkEquals(c(365, 1, 0), dim(ci@quantiles$tmax$inbase$q10))
  checkEquals(c(365, 2, 1), dim(ci2@quantiles$tmax$inbase$q10))

  ## Specifically test the adding and masking of NA values with threshold-based indices.
  checkEquals(structure(NA_real_, .Names = "1972-06"), climdex.tx10p(ci3)[30])
  
  tx90p.ci.valid <- structure(c(0, 0, 0, NA, NA, 0, 0, 0, 0, 0, 0, 0),
                              .Names = c("1970-01", "1970-02", "1970-03", "1970-04", "1970-05", "1970-06", "1970-07", "1970-08", "1970-09", "1970-10", "1970-11", "1970-12"))
  
  tx90p.ci2.valid <- structure(c(0, 0, 0, NA, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, 0, 0, 0, 0, 0, 0, 0),
                               .Names = c("1970-01", "1970-02", "1970-03", "1970-04", "1970-05", "1970-06", "1970-07", "1970-08", "1970-09", "1970-10", "1970-11", "1970-12", "1971-01", "1971-02", "1971-03", "1971-04", "1971-05", "1971-06", "1971-07", "1971-08", "1971-09", "1971-10", "1971-11", "1971-12"))
  
  checkEquals(tx90p.ci.valid, climdex.tx90p(ci))
  checkEquals(tx90p.ci2.valid, climdex.tx90p(ci2))
  
}

climdex.pcic.test.climdex.gsl <- function() {
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year", "jday")]), format="%Y %j", cal="gregorian")

  ## Load the data in.
  ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP, ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION, tmax.dates, tmin.dates, prec.dates, base.range=c(1981, 1990), northern.hemisphere=FALSE)

  checkEquals(get('ec.1018935.gsl.sh'), climdex.gsl(ci))
}

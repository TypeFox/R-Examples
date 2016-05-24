setOldClass("zoo")

setGeneric(
  name = 'phenoPhase',
  def = function(x, ...)
    standardGeneric("phenoPhase")
)

setMethod(
  f = "phenoPhase",
  signature = "ts",
  definition = function(x, mon.range = c(1, 12), ...) {

    d1 <- data.frame(yr = floor(time(x)), mon = cycle(x), val =
    	as.numeric(x))
    mons <- mon.range[1]:mon.range[2]
    d2 <- d1[d1$mon %in% mons, ]
    yrs <- unique(d2$yr)
    yrs.ok <- table(d2$yr, is.na(d2$val))[, 1] == length(mons)

    # max month
    a1 <- aggregate(d2$val, list(d2$yr), which.max)
    max.time <- ifelse(yrs.ok, a1$x, NA)
    max.time <- unlist(max.time) + mons[1] - 1

    # fulcrum
    fulc <- function(d) {
      if (sum(!is.na(d[,2])) < 2) {
        return(NA)
      } else {
        x <- d[, 1]
        y <- d[, 2]
        low <- min(x)
        up <- max(x)
        fun1 <- approxfun(x, y, rule = 2)
        fopt <- function(z) abs(integrate(fun1, low, z, ...)$value -
        	0.5*integrate(fun1, low, up, ...)$value)
        optimize(fopt, lower = low, upper = up)$minimum
      }
    }
    b1 <- by(d2[, c('mon', 'val')], as.factor(d2$yr), fulc)
    fulcrum <- round(ifelse(yrs.ok, as.numeric(b1), NA), 2)

    # weighted mean month
    weighted.mean.df <- function(d) {
      d <- na.omit(d)
      if (nrow(d) == 0) {
        return(NA)
      } else {
        weighted.mean(d[,1], d[,2])
      }
    }
    b2 <- by(d2[, c('mon', 'val')], as.factor(d2$yr), weighted.mean.df)
    mean.wt <- round(ifelse(yrs.ok, as.numeric(b2), NA), 2)

    as.data.frame(cbind(year = yrs, max.time, fulcrum, mean.wt),
    	row.names = 1:length(yrs))
  }
)

setMethod(
  f = "phenoPhase",
  signature = "zoo",
  definition = function(x, mon.range = c(1, 12), out = c('date', 'doy',
  	'julian'), ...) {

    # validate args
    if (!is(index(x), "DateTime"))
      stop('time index must be a DateTime object')
    indexx <- as.Date(index(x))
    out <- match.arg(out)

    d1 <- data.frame(julday = julian(indexx, origin =
    	as.Date("1970-01-01")), yr = years(indexx), mon =
    	monthNum(indexx), val = as.numeric(x))
    mons <- mon.range[1]:mon.range[2]
    d2 <- d1[d1$mon %in% mons, ]
    yrs <- unique(d2$yr)
    n <- table(d2$yr, is.na(d2$val))[, 1]

    # max day
    maxDay <- function(d)
      d[, 1][which.max(d[, 2])]
    b1 <- by(d2[, c('julday', 'val')], as.factor(d2$yr), maxDay)
    max.time <- as.numeric(b1)

    # fulcrum
    fulc <- function(d, m1 = mon.range[1], m2 = mon.range[2]) {
      if (sum(!is.na(d[,2])) < 2) {
        return(NA)
      } else {
        x <- d[, 1]
        y <- d[, 2]
        yr <- d[1, 3]
        mon.length <- c(31, ifelse(leapYear(yr), 29, 28), 31, 30, 31,
        	30, 31, 31, 30, 31, 30, 31)
        lo <- julian(as.Date(paste(yr, m1, 1, sep = '-')), origin =
        	as.Date("1970-01-01"))
        up <- julian(as.Date(paste(yr, m2, mon.length[m2], sep = '-')),
        	origin = as.Date("1970-01-01"))
        fun1 <- approxfun(x, y, rule = 2)
        fopt <- function(z) abs(integrate(fun1, lo, z, ...)$value -
        	0.5*integrate(fun1, lo, up, ...)$value)
        optimize(fopt, lower = lo, upper = up)$minimum
      }
    }
    b2 <- by(d2[, c('julday', 'val', 'yr')], as.factor(d2$yr), fulc)
    fulcrum <- ceiling(as.numeric(b2))

    # weighted mean
    weighted.mean.df <- function(d) {
      d <- na.omit(d)
      if (nrow(d) == 0) {
        return(NA)
      } else {
        weighted.mean(d[,1], d[,2])
      }
    }
    b3 <- by(d2[, c('julday', 'val')], as.factor(d2$yr),
    	weighted.mean.df)
    mean.wt <- ceiling(as.numeric(b3))

    switch(match(out, c('date', 'doy', 'julian')),
      data.frame(year = yrs, max.time = as.Date(max.time), fulcrum =
      	as.Date(fulcrum), mean.wt = as.Date(mean.wt), n, row.names =
      	NULL),
      data.frame(year = yrs, max.time =
      	as.POSIXlt(as.Date(max.time))$yday + 1, fulcrum =
      	as.POSIXlt(as.Date(fulcrum))$yday + 1, mean.wt =
      	as.POSIXlt(as.Date(mean.wt))$yday + 1, n, row.names = NULL),
      as.data.frame(cbind(year = yrs, max.time, fulcrum, mean.wt, n),
      	row.names = 1:length(yrs))
    )
  }
)

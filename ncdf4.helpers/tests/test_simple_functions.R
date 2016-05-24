library(PCICt)

ncdf4.helpers.test.nc.make.time.bounds <- function() {
  ts <- as.PCICt(c("1961-01-15", "1961-02-15", "1961-03-15"), cal="360")
  valid.data <- structure(c(-279936000, -277430400, -277344000, -274838400, -274752000,
                            -272246400, -272160000, -269654400, -269568000, -267062400, -266976000,
                            -264470400, -264384000, -261878400, -261792000, -259286400, -259200000,
                            -256694400, -256608000, -254102400, -254016000, -251510400, -251424000,
                            -248918400), .Dim = c(2L, 12L), cal = "360",
                          months = c(30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30), class = "PCICt", dpy = 360, tzone = "GMT", units = "secs")
  checkEquals(nc.make.time.bounds(ts, unit="month"), valid.data)
  checkException(nc.make.time.bounds(ts, unit="day"))
}

ncdf4.helpers.test.get.f.step.size <- function() {
  dat <- c(1, 2, 3, 4, 5, 7)
  checkEquals(get.f.step.size(dat, max), 2)
  checkEquals(get.f.step.size(dat, min), 1)
}

ncdf4.helpers.test.nc.is.regular.dimension <- function() {
  dat <- c(1, 2, 3, 4, 5, 6, 7)
  checkEquals(nc.is.regular.dimension(dat), TRUE)
  dat[7] <- 7.001
  checkEquals(nc.is.regular.dimension(dat), FALSE)
}

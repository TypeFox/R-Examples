library('smaa')

pvf1 <- function(x) {
  smaa.pvf(x, cutoffs=c(-0.15, 0.35), values=c(0, 1))
}
pvf2 <- function(x) {
  smaa.pvf(x, cutoffs=c(50, 100), values=c(1, 0))
}


stopifnot(all.equal(pvf1(0.35), 1))
stopifnot(all.equal(pvf1(-0.15), 0))
stopifnot(all.equal(pvf1(0.1), 0.5))
stopifnot(all.equal(pvf1(c(0.35, 0.1)), c(1, 0.5)))

stopifnot(all.equal(pvf2(50), 1))
stopifnot(all.equal(pvf2(100), 0))
stopifnot(all.equal(pvf2(75), 0.5))
stopifnot(all.equal(pvf2(c(75, 100)), c(0.5, 0)))

pvf3 <- function(x) {
  smaa.pvf(x, cutoffs=c(-0.15, 0.0, 0.25, 0.35), values=c(0, 0.1, 0.9, 1))
}
pvf4 <- function(x) {
  smaa.pvf(x, cutoffs=c(50, 75, 90, 100), values=c(1, 0.8, 0.5, 0))
}

stopifnot(all.equal(pvf3(0.35), 1.0))
stopifnot(all.equal(pvf3(-0.15), 0.0))
stopifnot(all.equal(pvf3(0.0), 0.1))
stopifnot(all.equal(pvf3(0.25), 0.9))
stopifnot(all.equal(pvf3(0.1), 2/5*0.8+0.1))

stopifnot(all.equal(pvf4(50), 1.0))
stopifnot(all.equal(pvf4(60), 1-(2/5*0.2)))
stopifnot(all.equal(pvf4(75), 0.8))
stopifnot(all.equal(pvf4(90), 0.5))
stopifnot(all.equal(pvf4(100), 0.0))
stopifnot(all.equal(pvf4(c(50,90,100)), c(1.0, 0.5, 0.0)))

pvf5 <- function(x) {
  smaa.pvf(x, cutoffs=c(50, 75, 90, 100), values=c(1, 0.8, 0.5, 0), outOfBounds="interpolate")
}

stopifnot(all.equal(pvf5(c(50,90,100)), c(1.0, 0.5, 0.0)))
stopifnot(all.equal(pvf5(c(10, 40, 101, 120)), c(1.32, 1.08, -0.05, -1)))

pvf6 <- function(x) {
  smaa.pvf(x, cutoffs=c(50, 75, 90, 100), values=c(1, 0.8, 0.5, 0), outOfBounds="clip")
}

stopifnot(all.equal(pvf6(c(50,90,100)), c(1.0, 0.5, 0.0)))
stopifnot(all.equal(pvf6(c(10, 40, 101, 120)), c(1.0, 1.0, 0.0, 0.0)))

pvf7 <- function(x) {
  smaa.pvf(x, cutoffs=c(-0.15, 0.0, 0.25, 0.35), values=c(0, 0.1, 0.9, 1), outOfBounds="clip")
}

stopifnot(all.equal(pvf7(c(0.35, -0.15, 0.0, 0.25)), c(1.0, 0.0, 0.1, 0.9)))
stopifnot(all.equal(pvf7(c(0.38, -0.20)), c(1.0, 0.0)))

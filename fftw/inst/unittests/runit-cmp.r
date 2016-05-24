##
## runit-cmp.r - Compare fftw and R core results
##

comparisonTest <- function(n) {
  ## x <- seq(-pi, pi, length.out=n)
  ## z <- sin(x) + cos(5*x)
  z <- runif(n)
  r <- fft(z)
  w <- FFT(z)
  checkEqualsNumeric(r, w)
  checkEqualsNumeric(z, Re(IFFT(w, scale=TRUE)))
}

test.3 <- function()
  comparisonTest(3)

test.4 <- function()
  comparisonTest(4)

test.5 <- function()
  comparisonTest(5)

test.128 <- function()
  comparisonTest(128)

test.2047 <- function()
  comparisonTest(2047)

test.2048 <- function()
  comparisonTest(2048)

test.2049 <- function()
  comparisonTest(2049)

test.262143 <- function()
    comparisonTest(262143)

test.262144 <- function()
    comparisonTest(262144)

test.262145 <- function()
    comparisonTest(262145)

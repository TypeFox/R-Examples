require("truncnorm")

################################################################################
## Check d/e/vtruncnorm all in one function:
check_dev <- function(a, b, mean=0, sd=1) {
  e <- etruncnorm(a, b, mean, sd)
  v <- vtruncnorm(a, b, mean, sd)

  ee <- integrate(function(x) x * dtruncnorm(x, a, b, mean, sd), a, b)$value
  if (abs(e - ee)/e > 0.00005) {
    message(sprintf("FAIL: etruncnorm(%4.1f, %4.1f, %4.1f, %4.1f) - mismatch %f vs. %f",
                    a, b, mean, sd, ee, e))
  }
  ev <- integrate(function(x) (x-ee)^2 * dtruncnorm(x, a, b, mean, sd), a, b)$value
  if (abs(v - ev)/v > 0.00005) {
    message(sprintf("FAIL: vtruncnorm(%4.1f, %4.1f, %4.1f, %4.1f) - mismatch %f vs. %f",
                    a, b, mean, sd, ev, v))
  }
}

## Left truncated:
check_dev(-3, Inf, 0, 1)
check_dev(-2, Inf, 1, 1)
check_dev( 2, Inf, 0, 1)
check_dev( 3, Inf, 1, 1)

check_dev(-3, Inf, 0, 2)
check_dev(-2, Inf, 1, 2)
check_dev( 2, Inf, 0, 2)
check_dev( 3, Inf, 1, 2)

## Doubly truncated:
check_dev(-3.0, -2.5, 0, 1)
check_dev(-3.0, -1.5, 0, 1)
check_dev(-3.0, -0.5, 0, 1)
check_dev(-3.0,  0.5, 0, 1)

check_dev(0.0, 0.5, 0, 1)
check_dev(0.0, 1.5, 0, 1)
check_dev(0.0, 2.5, 0, 1)
check_dev(0.0, 3.5, 0, 1)

## Extreme cases:
check_dev( 0.0, 1.0,  0.0, 10.0)
check_dev( 0.0, 1.0,  5.0,  1.0)
check_dev(-1.0, 0.0,  0.0, 10.0)
check_dev( 0.0, 1.0, -5.0,  1.0)

################################################################################
## Sanity checks on random number generators
check_r <- function(a, b, mean, sd, n=10000) {
  x <- rtruncnorm(n, a, b, mean, sd)
  if (!all(x > a & x < b)) {
    message(sprintf("FAIL: rtruncnorm(%i, %4.1f, %4.1f, %4.1f, %4.1f) - bounds",
                    n, a, b, mean, sd))
  }

  ## Check to make sure mean and variance have the correct magnitude.
  e.x <- mean(x)
  e <- etruncnorm(a, b, mean, sd)
  if (abs(e.x - e)/sd > 0.05) {
    message(sprintf("FAIL: rtruncnorm(%i, %4.1f, %4.1f, %4.1f, %4.1f) - mean %f vs. %f",
                    n, a, b, mean, sd, e.x, e))
  }
  sd.x <- sd(x)
  sd <- sqrt(vtruncnorm(a, b, mean, sd))
  if (abs(sd.x - sd)/sd.x > 0.05) {
    message(sprintf("FAIL: rtruncnorm(%i, %4.1f, %4.1f, %4.1f, %4.1f) - variance %f vs. %f",
                    n, a, b, mean, sd, sd.x, sd))
  }
}

## rtruncnorm == rnorm:
check_r(-Inf, Inf, 0, 1)

## 0 in (a, b):
check_r(-1, 1, 0, 1)
check_r(-1, 1, 1, 1)
check_r(-1, 1, 0, 2)

## 0 < (a, b):
check_r(1, 2, 0, 1)
check_r(1, 2, 1, 1)
check_r(1, 2, 0, 2)

## 0 > (a, b):
check_r(-2, -1, 0, 1)
check_r(-2, -1, 1, 1)
check_r(-2, -1, 0, 2)

## left truncation:
check_r(-2, Inf, 0, 1)
check_r(-2, Inf, 1, 1)
check_r(-2, Inf, 0, 2)
check_r( 0, Inf, 0, 1)
check_r( 0, Inf, 1, 1)
check_r( 0, Inf, 0, 2)
check_r( 2, Inf, 0, 1)
check_r( 2, Inf, 1, 1)
check_r( 2, Inf, 0, 2)

check_r(-0.2, Inf, 0, 1)
check_r(-0.2, Inf, 1, 1)
check_r(-0.2, Inf, 0, 2)
check_r( 0.0, Inf, 0, 1)
check_r( 0.0, Inf, 1, 1)
check_r( 0.0, Inf, 0, 2)
check_r( 0.2, Inf, 0, 1)
check_r( 0.2, Inf, 1, 1)
check_r( 0.2, Inf, 0, 2)

## Right truncation:
check_r(-Inf, -2, 0, 1)
check_r(-Inf, -2, 1, 1)
check_r(-Inf, -2, 0, 2)
check_r(-Inf,  0, 0, 1)
check_r(-Inf,  0, 1, 1)
check_r(-Inf,  0, 0, 2)
check_r(-Inf,  2, 0, 1)
check_r(-Inf,  2, 1, 1)
check_r(-Inf,  2, 0, 2)

check_r(-Inf, -0.2, 0, 1)
check_r(-Inf, -0.2, 1, 1)
check_r(-Inf, -0.2, 0, 2)
check_r(-Inf,  0.0, 0, 1)
check_r(-Inf,  0.0, 1, 1)
check_r(-Inf,  0.0, 0, 2)
check_r(-Inf,  0.2, 0, 1)
check_r(-Inf,  0.2, 1, 1)
check_r(-Inf,  0.2, 0, 2)

## Extreme examples:
check_r(-5, -4, 0, 1)

check_pq <- function(a, b, mean, sd) {
  for (p in runif(500)) {
    q <- qtruncnorm(p, a, b, mean, sd)
    pp <- ptruncnorm(q, a, b, mean, sd)
    if (abs(p - pp) > 0.00001) {
      message(sprintf("ptruncnorm(%6.4f, %6.4f, %6.4f, %6.4f, %6.4f) - disagree with qtruncnorm by %f",
                      p, a, b, mean, sd, abs(p - pp)))
    }
  }
}

check_pq(-1, 0, 0, 1)
check_pq(-1, 1, 0, 1)
check_pq( 1, 2, 0, 1)
check_pq(-1, 0, 4, 1)
check_pq(-1, 1, 4, 1)
check_pq( 1, 2, 4, 1)
check_pq(-1, 0, 0, 3)
check_pq(-1, 1, 0, 3)
check_pq( 1, 2, 0, 3)
check_pq(-1, Inf, 0, 1)
check_pq(-1, Inf, 4, 1)
check_pq(-1, Inf, 0, 3)
check_pq(-Inf, 1, 0, 1)
check_pq(-Inf, 1, 4, 1)
check_pq(-Inf, 1, 0, 3)

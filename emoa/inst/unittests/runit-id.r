##
## runit-pu.r - Pareto utility test
##

points <- matrix(c(1.0, 0.0, 0.0,
                   0.0, 1.0, 0.0,
                   0.0, 0.0, 1.0,
                   0.5, 0.5, 0.5,
                   0.5, 0.6, 0.6,
                   0.6, 0.5, 0.6,
                   0.6, 0.6, 0.5,
                   0.8, 0.8, 0.8),
                 ncol=8, byrow=FALSE)

nd <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
no <- c(1, 1, 1, 1, 2, 2, 2, 3)

test.hv_indicators <- function() {
  p1 <- points[, no==1]
  p2 <- points[, no==2]

  I12 <- hypervolume_indicator(p1, p2, ref=c(1, 1, 1))
  I21 <- hypervolume_indicator(p2, p1, ref=c(1, 1, 1))
  checkEquals(I21, 0.013)
  checkEquals(I12, -I21)

  I21p <- hypervolume_indicator(p2, p1, ref=c(10, 10, 10))
  checkTrue(I21 < I21p)
}

test.r_indicators <- function() {
  p1 <- points[, no==1]
  p2 <- points[, no==2]

  ## Basic sanity:
  checkEquals(r1_indicator(p1, p1), 0.5)
  checkEquals(r2_indicator(p1, p1), 0)
  checkEquals(r3_indicator(p1, p1), 0)  
  checkEquals(r1_indicator(p2, p2), 0.5)
  checkEquals(r2_indicator(p2, p2), 0)
  checkEquals(r3_indicator(p2, p2), 0)

  ## Precalculate values:
  r112 <- r1_indicator(p1, p2)
  r121 <- r1_indicator(p2, p1)

  r212 <- r2_indicator(p1, p2)
  r221 <- r2_indicator(p2, p1)

  r312 <- r3_indicator(p1, p2)
  r321 <- r3_indicator(p2, p1)
  
  ## Symmetry properties:
  checkEquals(r112 + r121, 1)
  checkEquals(r212, -r221)

  ## Known 'better':
  checkTrue(r112 > r121)
  checkTrue(r212 < r221)
  checkTrue(r312 < r321)
}

test.eps_indicator <- function() {
  p1 <- points[, no==1]
  p2 <- points[, no==2]

  checkEquals(epsilon_indicator(p1, p2), 0)
  checkEquals(epsilon_indicator(p2, p1), 0.6)
}

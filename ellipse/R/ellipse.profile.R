"ellipse.profile" <-
  function (x, which = c(1, 2), level = 0.95, t = sqrt(qchisq(level, 
                                                2)), npoints = 100, ...) 
{
  aa <- x[[which[1]]][[2]][, which[1]]
  ar <- x[[which[1]]][[2]][, which[2]]
  ra <- x[[which[2]]][[2]][, which[1]]
  rr <- x[[which[2]]][[2]][, which[2]]
  atau <- x[[which[1]]][[1]]
  rtau <- x[[which[2]]][[1]]
  arange <- range(c(aa, ra))
  rrange <- range(c(ar, rr))
  atau <- atau/t
  rtau <- rtau/t
  getad <- function(tau1, tau2) {
    if (abs(tau1) > 1) 
      tau1 <- tau1/abs(tau1)
    if (abs(tau2) > 1) 
      tau2 <- tau2/abs(tau2)
    acos1 <- acos(tau1)
    acos2 <- acos(tau2)
    d <- abs(acos1 - acos2)
    a <- (acos1 + acos2)/2
    if (acos1 < acos2) 
      a <- -a
    c(a, d)
  }
  myapprox <- function(x, y, where) {
    good <- is.finite(x) & is.finite(y)
    x <- x[good]
    y <- y[good]
    if (length(x) > 1) {
      result <- approx(x[good], y[good], where)$y
      bad <- is.na(result)
      if (any(bad)) {
        for (i in 1:length(result)) {
          if (bad[i]) {
            if (where[i] > x[length(x)]) {
              x1 <- x[length(x) - 1]
              y1 <- y[length(x) - 1]
              x2 <- x[length(x)]
              y2 <- y[length(x)]
            }
            else if (where[i] < x[1]) {
              x1 <- x[1]
              y1 <- y[1]
              x2 <- x[2]
              y2 <- y[2]
            }
            else stop("Unexpected NA")
            result[i] <- y1 + (where[i] - x1)/(x2 - x1) * 
              (y2 - y1)
          }
        }
      }
    }
    else result <- rep(y, length(where))
    result
  }
  ad <- matrix(NA, nrow = 5, ncol = 2)
  ad[1, ] <- getad(1, myapprox(rr, rtau, myapprox(aa, ar, myapprox(atau, aa, 1))))
  ad[2, ] <- getad(myapprox(aa, atau, myapprox(rr, ra, myapprox(rtau, rr, 1))), 1)
  ad[3, ] <- getad(-1, myapprox(rr, rtau, myapprox(aa, ar, myapprox(atau, aa, -1))))
  ad[4, ] <- getad(myapprox(aa, atau, myapprox(rr, ra, myapprox(rtau, rr, -1))), -1)
  i <- order(ad[1:4, 1])
  ad[1:4, ] <- ad[i, ]
  ad[5, 1] <- ad[1, 1] + 2 * pi
  ad[5, 2] <- ad[1, 2]
  ad <- ad[!duplicated(ad[, 1]), ]
  adv <- spline(ad, n = npoints, method= "periodic")
  avals <- adv$x
  dvals <- adv$y
  matrix(c(myapprox(atau, aa, cos(avals + dvals/2)), myapprox(rtau,  rr, cos(avals - dvals/2))), length(avals), 2, dimnames = list(NULL, names(x[which])))
}

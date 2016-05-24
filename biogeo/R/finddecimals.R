finddecimals <-
function (dat, x = "x", y = "y") 
{
  cn <- names(dat)
  f1 <- match(x, cn)
  f2 <- match(y, cn)
  nr <- nrow(dat)
  x <- as.character(dat[, f1])
  y <- as.character(dat[, f2])
  b <- rep(5, nr)
  xv <- b
  for (i in 1:nr) {
    suppressWarnings(a <- as.numeric(x[i]))
    if (is.na(a)) {
      b[i] <- 0
    }
    else {
      a1 <- abs(a) #- floor(abs(a))
      b[i] <- ifelse(a1 > 0, 1, 0)
      xv[i] <- ifelse(a1 > 0, a, 0)
    }
  }
  bx <- b
  b <- rep(5, nr)
  yv <- b
  for (i in 1:nr) {
    suppressWarnings(a <- as.numeric(y[i]))
    if (is.na(a)) {
      b[i] <- 0
    }
    else {
      a1 <- abs(a) #- floor(abs(a))
      b[i] <- ifelse(a1 > 0, 1, 0)
      yv[i] <- ifelse(a1 > 0, a, 0)
    }
  }
  by <- b
  x1 <- (abs(xv) > 180) * 1
  y1 <- (abs(yv) > 90) * 1
  fx1 <- which(x1 == 1)
  fy1 <- which(y1 == 1)
  bx[fx1] <- 0
  by[fy1] <- 0
  decs <- bx * by
  return(decs)
}

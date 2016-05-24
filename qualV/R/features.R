standard01 <- function (x) {
  xn <- (x - min(x)) / (max(x) - min(x))
  xn
}

steps <- function (x, difx, tol = 1e-6 * max(abs(difx))) {
  if (length(x) > 2 && any(abs(diff(difx)) > tol))
    warning ("time steps may be unequal", immediate. = TRUE)
}

scaledEpsilon <- function(difq, scale) {
  epsilon <- switch(scale,
    mean   = mean(abs(difq[-1])),
    range  = abs(diff(range(difq))),
    IQR    = abs(diff(IQR(difq))),
    sd     = sd(difq[-1]),
    none   = 1
  )
}

f.slope <- function (x, y, f = 0.1, scale = c("mean", "range", "IQR", "sd", "none")) {
  N(x, y)
  difx <- diff(x)
  steps(x, difx) # checks for equal spacing
  dify <- diff(y)
  difq <- c(0, dify / difx)

  ## rescale differences
  scale   <- match.arg(scale)
  epsilon <- scaledEpsilon(difq, scale) * f
  
  dat <- data.frame(difq, v = rep(0, length(difq)))
  ## increase
  plus <- which(dat$difq > epsilon)
  dat$v[plus] <- "A"
  ## decrease
  minus <- which(dat$difq < -epsilon)
  dat$v[minus] <- "B"
  ## constant
  null <- which(dat$difq >= -epsilon & dat$difq <= epsilon)

  dat$v[null] <- "C"
  dat$v
}

f.curve <- function (x, y, f = 0.1, scale = c("mean", "range", "IQR", "sd", "none")) {
  N(x, y)
  
  difx1 <- diff(x)
  steps(x, difx1)
  
  dify1 <- diff(y)
  difq1 <- dify1 / difx1

  difx2 <- diff(x, lag = 2)
  dify2 <- diff(difq1)
  difq2 <- dify2 / difx2
  
  dat <- data.frame(difq2 = difq2, v = rep(0, length(difq2)))

  ## rescale differences
  scale   <- match.arg(scale)
  epsilon <- scaledEpsilon(difq2, scale) * f
  
  # convex
  plus <- which(dat$difq2 > epsilon)
  dat$v[plus] <- "K"
  # concave
  minus <- which(dat$difq2 < -epsilon)
  dat$v[minus] <- "I"
  # constant
  null <- which(dat$difq2 >= -epsilon & dat$difq2 <= epsilon)
  dat$v[null] <- "J"
  
  v <- c(0, paste(c(0, dat$v), c(dat$v, 0), sep = ""))
  v
}

f.steep <- function (x, y, f1 = 1, f2 = 0.1) {
  N(x, y)
  
  y <- standard01(y)
  difx <- diff(x)
  steps(x, difx)
  
  dify <- diff(y)
  difq <- c(0, dify / difx)
  alpha <- abs(atan(difq)) * 180 / pi
  dat <- data.frame(difq = difq, alpha = alpha, v = rep(0, length(difq)))
     
  # very steep
  ss <- which(dat$alpha > f1)
  dat$v[ss] <- "S"
  # steep
  steep <- which(dat$alpha >= f2 & dat$alpha <= f1)
  dat$v[steep] <- "T"
  # not steep
  nsteep <- which(dat$alpha < f2)
  dat$v[nsteep] <- "U"
  
  v <- dat$v
  v
}

f.level <- function (y, high = 0.8, low = 0.2) {
  if(length(y) == 0) stop ("vector of length zero")

  y <- standard01(y)
  dat <- data.frame(y, v = rep(0, length(y)))
  
  H <- which(dat$y >= high)
  dat$v[H] <- "H"
  M <- which(dat$y > low & dat$y < high)
  dat$v[M] <- "M"
  N <- which(dat$y <= low)
  dat$v[N] <- "L"

  v <- dat$v
  v
}

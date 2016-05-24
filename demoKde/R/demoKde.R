kde <- function(x, bw = bw.nrd0, kernel = kernelGaussian, n = 4096,
                from = min(x) - 3*sd, to = max(x) + 3*sd, adjust = 1,
                ...) {
  if(has.na <- any(is.na(x))) {
    x <- na.omit(x)
    if(length(x) == 0)
        stop("no finite or non-missing data!")
  }
  sd <- (if(is.numeric(bw)) bw[1] else bw(x)) * adjust
  X <- seq(from, to, len = n)
  M <- outer(X, x, kernel, sd = sd, ...)
  structure(list(x = X, y = rowMeans(M), bw = sd,
                 call = match.call(), n = length(x),
                 data.name = deparse(substitute(x)),
                 has.na = has.na), class =  "density")
}

kernelBiweight <- function(x, mean=0, sd=1) {
  h <- sqrt(7)*sd
  ifelse((z <- abs(x-mean)) < h, 15/16*(1 - (z/h)^2)^2/h, 0)
}

kernelCosine <- function(x, mean=0, sd=1) {
  h <- sqrt(1/(1-8/pi^2))*sd
  ifelse((z <- abs(x-mean)) < h, pi/4*cos((pi*z)/(2*h))/h, 0)
}

kernelEpanechnikov <- function(x, mean=0, sd=1) {
  h <- sqrt(5)*sd
  ifelse((z <- abs(x-mean)) < h, 3/4*(1 - (z/h)^2)/h, 0)
}

kernelGaussian <- function(x, mean=0, sd=1)
    dnorm(x, mean = mean, sd = sd)

kernelLogistic <- function(x, mean=0, sd=1)
    stats::dlogis(x, mean, sqrt(3)/pi*sd)

kernelOptCosine <- function(x, mean=0, sd=1) {
  h <- sqrt(1/(1-8/pi^2))*sd
  ifelse((z <- abs(x-mean)) < h, pi/4*cos((pi*z)/(2*h))/h, 0)
}

kernelRectangular <- function(x, mean=0, sd=1) {
  h <- sqrt(3)*sd
  ifelse(abs(x-mean) < h, 1/(2*h), 0)
}

kernelSquaredCosine <- function(x, mean=0, sd=1) {
  h <- sqrt(3/(1-6/pi^2))*sd
  ifelse((z <- abs(x-mean)) < h, cos(pi*z/(2*h))^2/h, 0)
}

kernelTriangular <- function(x, mean=0, sd=1) {
  h <- sqrt(24)*sd/2
  ifelse((z <- abs(x-mean)) < h, (1 - z/h)/h, 0)
}

kernelTricube <- function(x, mean=0, sd=1) {
  h <- sqrt(243/35)*sd
  ifelse((z <- abs(x - mean)) < h, 70/81*(1 - (z/h)^3)^3/h, 0)
}

kernelTriweight <- function(x, mean=0, sd=1) {
  h <- sqrt(9)*sd
  ifelse((z <- abs(x-mean)) < h, 35/32*(1 - (z/h)^2)^3/h, 0)
}

kernelUniform <- function(x, mean=0, sd=1) {
  h <- sqrt(3)*sd
  ifelse(abs(x-mean) < h, 1/(2*h), 0)
}

"h.select" <- function(x, y = NA, weights = NA, group = NA, ...) {

  data    <- sm.check.data(x, y, weights = weights, group = group, ...)
  x       <- data$x
  y       <- data$y
  weights <- data$weights
  group   <- data$group
  nobs    <- data$nobs
  ndim    <- data$ndim
  density <- data$density
  opt     <- data$options
  
  if (all(!is.na(group))) {
     group.fac <- factor(group)
     h.all <- matrix(0, ncol = ndim, nrow = 0)
     for (igroup in 1:length(levels(group.fac))) {
        level.i <- levels(group.fac)[igroup]
        if (ndim == 1)
           h.igroup <- h.select(x[group.fac == level.i],
                                y[group.fac == level.i],
                                weights[group.fac == level.i], ...)
        else
           h.igroup <- h.select(x[group.fac == level.i,],
                                y[group.fac == level.i],
                                weights[group.fac == level.i], ...)
        h.all <- rbind(h.all, h.igroup)
        }
    h.gmean <- apply(h.all, 2, FUN = function(x) exp(mean(log(x))))
    return(as.vector(h.gmean))
    }

  if (ndim == 1) 
     replace.na(opt, df, 6)
  else if (ndim == 2)
     replace.na(opt, df, 12)
  if ((!density) & opt$df <= 2)
     stop("df must be > 2")
  if (density)
     replace.na(opt, method, "normal")
  else
     replace.na(opt, method, "df")
  if ((ndim == 3) & !(density & opt$method == "normal"))
     stop("bandwidth selection not available for 3 dimensions.")

  method       <- opt$method
  df           <- opt$df
  structure.2d <- opt$structure.2d
  if (method == "df" & ndim == 2 & structure.2d == "different")
     stop("df method is not appropriate for different bandwidths")

  replace.na(opt, nbins, round((nobs > 100) * 8 * log(nobs) / ndim))
  if (opt$nbins > 0 & ndim < 3) {
     if (!all(weights == 1) & opt$verbose > 0) 
            cat("Warning: weights overwritten by binning\n")
     data <- binning(x, y, nbins = opt$nbins)
     }
  else
     data <- list(x = x, means = y, x.freq = rep(1,nobs),
                  devs = rep(0,nobs))

  h.weights <- opt$h.weights
  if (opt$verbose > 0 & !all(is.na(h.weights)) & opt$nbins > 0) {
     h.weights <- rep(1, length(data$x.freq))
     cat("h.weights cannot be used with binning")
     }
  if (all(is.na(h.weights))) h.weights <- rep(1, length(data$x.freq))
  data$h.weights <- h.weights

  sd <- sqrt(diag(as.matrix(var(x))))
  if (ndim==1)
     start <- sd / 2
  else
     start <- switch(structure.2d, common   = mean(sd / 2),
                                   scaled   = 0.5, 
                                   separate = sd / 2)
  data$sd <- sd

  if (density & method == "normal")
      return(hnorm(x))
  else if (density & method == "sj") {
      if (ndim > 1) stop("Sheather-Jones method requires 1-d data")
      if (!all(weights == 1) & opt$verbose > 0)
            cat("Warning: weights are not used in the sj method\n")
      return(hsj(x))
      }
  else {
     if (density) crit.type <- "dens" else crit.type <- "reg"
     fname <- paste(method, ".crit", ".", crit.type, sep = "")
     if (structure.2d == "separate") {
        result <- optim(par = log(start), fn = get(fname),
                        control = list(reltol = 1e-6),
                        data = data, structure.2d = structure.2d, opt = opt)
        h.result <- exp(result$par)
        }
     else {
        result <- optimise(f = get(fname),
                           interval = log(c(start / 8, start * 4)),
                           data = data, structure.2d = structure.2d, opt = opt)
        h.result <- exp(result$minimum)
        }
     }

  if (ndim == 2)
    h.result <- switch(structure.2d, common = rep(h.result,2),
                                     scaled = h.result * sd,
                                     h.result)

  return(h.result)

  }

"cv" <- function (x, h, ...) {
    opt <- sm.options(list(...))
    if (!isMatrix(x)) {
        n <- length(x)
        replace.na(opt, h.weights, rep(1, n))
        hcvff <- sum(dnorm(0, mean = 0, sd = sqrt(2) * h * opt$h.weights))/(n *
            (n - 1))
        W <- matrix(rep(x, rep(n, n)), ncol = n, byrow = TRUE)
        W <- W - matrix(rep(x, n), ncol = n, byrow = TRUE)
        W1 <- matrix(rep(opt$h.weights^2, n), ncol = n, byrow = TRUE)
        W2 <- exp(-0.5 * (W/(h * sqrt(W1 + t(W1))))^2)/(sqrt(2 * pi) *
                                                        h * sqrt(W1 + t(W1)))
        hcvff <- hcvff + (sum(W2) - sum(diag(W2))) * (n - 2)/(n * (n - 1)^2)
        W2 <- exp(-0.5 * (W/(h * sqrt(W1)))^2)/(sqrt(2 * pi) * h * sqrt(W1))
        hcvff <- hcvff - (sum(W2) - sum(diag(W2))) * 2/(n * (n - 1))
        }
    if (isMatrix(x)) {
        x1 <- x[, 1]
        x2 <- x[, 2]
        h1 <- h * sqrt(var(x1))
        h2 <- h * sqrt(var(x2))
        n <- length(x1)
        replace.na(opt, h.weights, rep(1, n))
        hcvff <- sum(dnorm(0, mean = 0, sd = sqrt(2) * h1 * opt$h.weights) *
            dnorm(0, mean = 0, sd = sqrt(2) * h2 * opt$h.weights))/(n *
            (n - 1))
        W <- matrix(rep(x1, rep(n, n)), ncol = n, byrow = TRUE)
        W <- W - matrix(rep(x1, n), ncol = n, byrow = TRUE)
        W1 <- matrix(rep(opt$h.weights^2, n), ncol = n, byrow = TRUE)
        W2 <- exp(-0.5 * (W/(h1 * sqrt(W1 + t(W1))))^2)/(sqrt(2 *
            pi) * h1 * sqrt(W1 + t(W1)))
        W <- matrix(rep(x2, rep(n, n)), ncol = n, byrow = TRUE)
        W <- W - matrix(rep(x2, n), ncol = n, byrow = TRUE)
        W2 <- W2 * exp(-0.5 * (W/(h2 * sqrt(W1 + t(W1))))^2)/(sqrt(2 *
            pi) * h2 * sqrt(W1 + t(W1)))
        hcvff <- hcvff + (sum(W2) - sum(diag(W2))) * (n - 2)/(n *
            (n - 1)^2)
        W2 <- exp(-0.5 * (W/(h2 * sqrt(W1)))^2)/(sqrt(2 * pi) *
            h2 * sqrt(W1))
        W <- matrix(rep(x1, rep(n, n)), ncol = n, byrow = TRUE)
        W <- W - matrix(rep(x1, n), ncol = n, byrow = TRUE)
        W2 <- W2 * exp(-0.5 * (W/(h1 * sqrt(W1)))^2)/(sqrt(2 *
            pi) * h1 * sqrt(W1))
        hcvff <- hcvff - (sum(W2) - sum(diag(W2))) * 2/(n * (n - 1))
        }
    hcvff
    }

"hcv" <- function (x, y = NA, hstart = NA, hend = NA, ...) {
    opt <- sm.options(list(...))
    replace.na(opt, ngrid, 8)
    replace.na(opt, display, "none")
    if (length(dim(x)) > 0) {
        ndim <- 2
        n <- length(x[, 1])
        }
    else {
        ndim <- 1
        n <- length(x)
        }
    replace.na(opt, h.weights, rep(1, n))
    ngrid <- opt$ngrid
    display <- opt$display
    h.weights <- opt$h.weights
    if (is.na(hstart)) {
        if (ndim == 1) hstart <- hnorm(x)/10
        else if (any(is.na(y)))
            hstart <- hnorm(x[, 1]/sqrt(var(x[, 1])))/10
        else hstart <- hnorm(x[, 1]/sqrt(var(x[, 1])))/4
        }
    if (is.na(hend)) {
        if (ndim == 1)
            hend <- hnorm(x) * 2
        else hend <- hnorm(x[, 1]/sqrt(var(x[, 1]))) * 2
        }
    cvgrid <- vector("numeric", length = ngrid)
    hgrid <- log(hstart) + (log(hend) - log(hstart)) * (0:(ngrid - 1)) / (ngrid - 1)
    if (any(is.na(y))) {
        for (i in 1:ngrid) cvgrid[i] <- cv(x, exp(hgrid[i]),
            h.weights = h.weights)
        }
    else {
        if (ndim == 1)
            for (i in 1:ngrid) {
                cvgrid[i] <- sum((y - sm.weight(x, x, h = exp(hgrid[i]),
                  cross = TRUE, options = list(h.weights = h.weights)) %*% y)^2)
                }
        if (ndim == 2)
            for (i in 1:ngrid) {
                cvgrid[i] <- sum((y - sm.weight2(x, x, 
                        exp(hgrid[i] * c(sqrt(var(x[, 1])), sqrt(var(x[, 2])))), 
                        cross = TRUE, options = list(h.weights = h.weights)) %*% y)^2)
                }
        }
    if (any(is.na(cvgrid))) {
        cat("\n")
        cat("hcv: some computations failed.", "\n")
        cat("Try readjusting hstart and hend.", "\n")
        cat("hstart: ", hstart, "\n")
        cat("hend  : ", hend, "\n")
        cat("\n")
        print(cbind(h = exp(hgrid), cv = cvgrid))
        stop()
        }
    ind <- (1:ngrid)[cvgrid == min(cvgrid)]
    if (!(display == "none")) {
        if (!opt$add) {
            if (display == "log")
                plot(hgrid, cvgrid, type = "l", xlab = "Log h", ylab = "CV")
            else plot(exp(hgrid), cvgrid, type = "l", xlab = "h", ylab = "CV")
            }
        else {
            if (display == "log")
                lines(hgrid, cvgrid)
            else lines(exp(hgrid), cvgrid)
            }
        }
    if (ind == 1 | ind == ngrid) {
        cat("\n")
        cat("hcv: boundary of search area reached.", "\n")
        cat("Try readjusting hstart and hend.", "\n")
        cat("hstart: ", hstart, "\n")
        cat("hend  : ", hend, "\n")
        cat("\n")
        print(cbind(h = exp(hgrid), cv = cvgrid))
        stop()
        }
    v0 <- cvgrid[ind - 1]
    v1 <- cvgrid[ind]
    v2 <- cvgrid[ind + 1]
    l0 <- hgrid[ind - 1]
    l1 <- hgrid[ind]
    l2 <- hgrid[ind + 1]
    aa <- (v1 - v0 - (l1 - l0) * (v1 - v2)/(l1 - l2))/(l1^2 -
        l0^2 - (l1^2 - l2^2) * (l1 - l0)/(l1 - l2))
    bb <- (v1 - v2 - aa * (l1^2 - l2^2))/(l1 - l2)
    cc <- v0 - aa * l0^2 - bb * l0
    h <- exp(-bb/(2 * aa))
    if (ndim == 1) result <- h
    else result <- c(h * sqrt(var(x[, 1])), h * sqrt(var(x[, 2])))
    result
    }

"hnorm" <- function (x, weights = NA) {
    if (isMatrix(x)) {
        if (all(is.na(weights)))
            weights <- rep(1, nrow(x))
        ndim <- ncol(x)
        n <- sum(weights)
        sd <- sqrt(apply(x, 2, wvar, w = weights))
        if (ndim == 2)
            hh <- sd * (1/n)^(1/6)
        else if (ndim == 3)
            hh <- sd * (4/(5 * n))^(1/7)
        if (ndim > 3) stop("data with >3 dimensions are not allowed.")
        hh
        }
    else {
        if (all(is.na(weights)))
            weights <- rep(1, length(x))
        sd <- sqrt(wvar(x, weights))
        sd * (4/(3 * sum(weights)))^(1/5)
        }
    }

"hsj" <- function (x) {
    h0 <- hnorm(x)
    v0 <- sj(x, h0)
    if (v0 > 0) hstep <- 1.1
    else hstep <- 0.9
    h1 <- h0 * hstep
    v1 <- sj(x, h1)
    while (v1 * v0 > 0) {
        h0 <- h1
        v0 <- v1
        h1 <- h0 * hstep
        v1 <- sj(x, h1)
        }
    h0 + (h1 - h0) * abs(v0)/(abs(v0) + abs(v1))
    }


"sj" <- function (x, h) {
    phi6 <- function(x) (x^6 - 15 * x^4 + 45 * x^2 - 15) * dnorm(x)
    phi4 <- function(x) (x^4 - 6 * x^2 + 3) * dnorm(x)
    n <- length(x)
    lambda <- quantile(x, 0.75) - quantile(x, 0.25)
    a <- 0.92 * lambda * n^(-1/7)
    b <- 0.912 * lambda * n^(-1/9)
    W <- matrix(rep(x, rep(n, n)), ncol = n, byrow = TRUE)
    W <- W - matrix(rep(x, n), ncol = n, byrow = TRUE)
    W1 <- matrix(phi6(W/b), ncol = n)
    tdb <- as.numeric(rep(1, n) %*% W1 %*% rep(1, n))
    tdb <- -tdb/(n * (n - 1) * b^7)
    W1 <- matrix(phi4(W/a), ncol = n)
    sda <- as.numeric(rep(1, n) %*% W1 %*% rep(1, n))
    sda <- sda/(n * (n - 1) * a^5)
    alpha2 <- 1.357 * (abs(sda/tdb))^(1/7) * h^(5/7)
    W1 <- matrix(phi4(W/alpha2), ncol = n)
    sdalpha2 <- as.numeric(rep(1, n) %*% W1 %*% rep(1, n))
    sdalpha2 <- sdalpha2/(n * (n - 1) * alpha2^5)
    result <- (dnorm(0, sd = sqrt(2))/(n * abs(sdalpha2)))^0.2 - h
    attributes(result)$names <- NULL
    as.double(result)
    }


"df.crit.reg" <- function(log.h, data, structure.2d, opt) {
  x    <- data$x
  freq <- data$x.freq
  h    <- exp(log.h)
  if (is.vector(x))
     S <- sm.weight(x, x, h, weights = freq, options = opt)
  if (is.matrix(x)) {
     h <- switch(structure.2d, scaled = h * data$sd, common = rep(h, 2), h)
     S <- sm.weight2(x, x, h, weights = freq, options = opt)
     }
  (sum(diag(S)) - opt$df)^2
  }

"cv.crit.reg" <- function(log.h, data, structure.2d, opt) {
   x    <- data$x
   y    <- data$means
   freq <- data$x.freq
   h    <- exp(log.h)
   if (is.vector(x))
      S <- sm.weight(x, x, h, cross = T, weights = freq, options = opt)
   if (is.matrix(x)) {
      h <- switch(structure.2d, scaled = h * data$sd, common = rep(h, 2), h)
      S <- sm.weight2(x, x, h, cross = T, weights = freq, options = opt)
      }
   n  <- length(y)
   cv <- sum(freq * ((diag(n) - S) %*% y)^2 + data$devs) / sum(freq)
   if(opt$verbose > 1) cat("h, CV: ", h, cv, "\n")
   if(is.na(cv)) cv <- 1e10
   cv
   }

"aicc.crit.reg" <- function(log.h, data, structure.2d, opt) {
   x    <- data$x
   y    <- data$means
   freq <- data$x.freq
   h    <- exp(log.h)
   if (is.vector(x))
     S <- sm.weight(x, x, h, weights=freq, options=opt)
   if (is.matrix(x)) {
     h <- switch(structure.2d, scaled = h * data$sd, common = rep(h, 2), h)
     S <- sm.weight2(x, x, h, weights=freq, options=opt)
     }
   tr.S <- sum(freq * diag(S))
   nobs <- sum(freq)
   n <- length(y)
   sig.sq <-  sum(freq * ((diag(n) - S) %*% y)^2 + data$devs) / nobs 
   penalty <- 1 + 2 * (tr.S + 1) / (nobs - tr.S - 2)
   aicc <- log(sig.sq) + penalty
   if(opt$verbose > 1) cat("h, AIC.c: ", h, aicc,"\n")
   if(is.na(aicc)) aicc <- 1e10
   aicc
   }

"cv.crit.dens" <- function(log.h, data, structure.2d, opt) {
   h         <- exp(log.h)
   x         <- data$x
   freq      <- data$x.freq
   n         <- length(freq)
   h.weights <- data$h.weights
   if (!is.matrix(x)) {
        hcvff <- sum(freq * dnorm(0, 0, sqrt(2)*h*h.weights))/(n*(n-1))
    W     <- matrix(rep(x, n), ncol = n)
    W     <- W - t(W)
    W1    <- matrix(rep(h.weights^2, n), ncol = n, byrow = TRUE)
    W2    <- exp(-.5 * (W/(h*sqrt(W1+t(W1))))^2) /
            (sqrt(2*pi)*h*sqrt(W1+t(W1)))
    W2    <- W2 * matrix(rep(freq,n), ncol=n) *
                  matrix(rep(freq,n), ncol=n, byrow = TRUE)
    hcvff <- hcvff + (sum(W2) - sum(diag(W2)))*(n-2)/(n*(n-1)^2)
    W2    <- exp(-.5 * (W/(h*sqrt(W1)))^2) / (sqrt(2*pi)*h*sqrt(W1))
    W2    <- W2 * matrix(rep(freq,n),ncol=n) *
                  matrix(rep(freq,n),ncol=n, byrow = TRUE)
    hcvff <- hcvff - (sum(W2) - sum(diag(W2)))*2/(n*(n-1))
    }

   if (is.matrix(x)) {
        h <- switch(structure.2d, scaled = h * data$sd,
                              common = rep(h, 2),    h)
    x1 <- x[,1]
    x2 <- x[,2]
    h1 <- h[1]
    h2 <- h[2]
    hcvff <- sum(freq * dnorm(0, 0, sqrt(2) * h1 * h.weights) *
            dnorm(0, 0, sqrt(2) * h2 * h.weights))/(n*(n-1))
    W     <- matrix(rep(x1, n), ncol = n)
    W     <- W - t(W)
    W1    <- matrix(rep(h.weights^2, n),  ncol = n, byrow = TRUE)
    W2    <- exp(-.5 * (W/(h1 * sqrt(W1+t(W1))))^2) /
            (sqrt(2 * pi) * h1 * sqrt(W1+t(W1)))
    W     <- matrix(rep(x2, n), ncol = n)
    W     <- W - t(W)
    W2    <- W2 * exp(-.5 * (W/(h2 * sqrt(W1+t(W1))))^2) /
            (sqrt(2 * pi) * h2 * sqrt(W1+t(W1)))
    W2    <- W2 * matrix(rep(freq,n), ncol=n) *
                  matrix(rep(freq,n), ncol=n, byrow = TRUE)
    hcvff <- hcvff + (sum(W2) - sum(diag(W2)))*(n-2)/(n*(n-1)^2)

    W2    <- exp(-.5 * (W/(h2 * sqrt(W1)))^2) /
            (sqrt(2 *pi) * h2 * sqrt(W1))
    W     <- matrix(rep(x1, n), ncol = n)
    W     <- W - t(W)
    W2    <- W2 * exp(-.5 * (W/(h1 * sqrt(W1)))^2) /
            (sqrt(2 * pi) * h1 * sqrt(W1))
    W2    <- W2 * matrix(rep(freq,n), ncol=n) *
                  matrix(rep(freq,n), ncol=n, byrow = TRUE)
    hcvff <- hcvff - (sum(W2) - sum(diag(W2))) * 2 / (n*(n-1))

    }
    hcvff
}

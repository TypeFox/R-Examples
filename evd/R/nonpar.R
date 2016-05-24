
"abvnonpar"<- 
function(x = 0.5, data, epmar = FALSE, nsloc1 = NULL, nsloc2 = NULL,
  method = c("cfg","pickands","tdo","pot"), k = nrow(data)/4,
  convex = FALSE, rev = FALSE, madj = 0, kmar = NULL, plot = FALSE,
  add = FALSE, lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1,
  xlim = c(0,1), ylim = c(0.5,1), xlab = "t", ylab = "A(t)", ...)
{
  if(mode(x) != "numeric" || any(x < 0, na.rm=TRUE) ||
    any(x > 1, na.rm=TRUE)) stop("invalid argument for `x'")
  method <- match.arg(method)

  # Empirical transform to exponential
  epdata <- apply(data, 2, rank, na.last = "keep")
  nasm <- apply(data, 2, function(x) sum(!is.na(x)))
  epdata <- epdata / rep(nasm+1, each = nrow(data))
  epdata <- -log(epdata)
  if(epmar) data <- epdata
  # End empirical transform  

  if(!epmar) {

    if(method == "pot") {
      # Parametric pot transform to exponential
      if(any(k >= nasm)) stop("k is too large")
      u1 <- sort(data[,1], decreasing = TRUE)[k+1]
      u2 <- sort(data[,2], decreasing = TRUE)[k+1]
      d1ab <- (data[,1] > u1) & !is.na(data[,1])
      d2ab <- (data[,2] > u2) & !is.na(data[,2])
      if(!is.null(kmar)) {
        data[d1ab,1] <- mtransform(data[d1ab,1], c(u1, kmar))
        data[d2ab,2] <- mtransform(data[d2ab,2], c(u2, kmar))
      }
      else {
        mle.m1 <- c(u1, fitted(fpot(data[d1ab,1], threshold = u1)))
        mle.m2 <- c(u2, fitted(fpot(data[d2ab,2], threshold = u2)))
        data[d1ab,1] <- mtransform(data[d1ab,1], mle.m1)
        data[d2ab,2] <- mtransform(data[d2ab,2], mle.m2)   
      }
      data[d1ab,1] <- -log(1 - k * data[d1ab,1] / nasm[1])
      data[d2ab,2] <- -log(1 - k * data[d2ab,2] / nasm[2])
      data[!d1ab, 1] <- epdata[!d1ab, 1]
      data[!d2ab, 2] <- epdata[!d2ab, 2]
      # End parametric pot transform
    }

    if(method != "pot") {
      # Parametric gev transform to exponential
      if(!is.null(kmar)) {
        data <- mtransform(data, kmar)
      }
      else {
        if(!is.null(nsloc1)) {
          if(is.vector(nsloc1)) nsloc1 <- data.frame(trend = nsloc1)
          if(nrow(nsloc1) != nrow(data))
            stop("`nsloc1' and data are not compatible")
          nslocmat1 <- cbind(1,as.matrix(nsloc1))
        }
        if(!is.null(nsloc2)) {
          if(is.vector(nsloc2)) nsloc2 <- data.frame(trend = nsloc2)
          if(nrow(nsloc2) != nrow(data))
            stop("`nsloc2' and data are not compatible")
          nslocmat2 <- cbind(1,as.matrix(nsloc2))
        }
        mle.m1 <- fitted(fgev(data[,1], nsloc = nsloc1, std.err = FALSE))
        loc.mle.m1 <- mle.m1[grep("^loc", names(mle.m1))]
        if(is.null(nsloc1)) loc.mle.m1 <- rep(loc.mle.m1, nrow(data))
        else loc.mle.m1 <- nslocmat1 %*% loc.mle.m1
        mle.m1 <- cbind(loc.mle.m1, mle.m1["scale"], mle.m1["shape"])
        mle.m2 <- fitted(fgev(data[,2], nsloc = nsloc2, std.err = FALSE))
        loc.mle.m2 <- mle.m2[grep("^loc", names(mle.m2))]
        if(is.null(nsloc2)) loc.mle.m2 <- rep(loc.mle.m2, nrow(data))
        else loc.mle.m2 <- nslocmat2 %*% loc.mle.m2
        mle.m2 <- cbind(loc.mle.m2, mle.m2["scale"], mle.m2["shape"])
        data <- mtransform(data, list(mle.m1, mle.m2))
      }
      # End parametric gev transform
    }
  }
    
  if(rev) data <- data[,2:1]
  data <- na.omit(data)
  if(plot || add) x <- seq(0, 1, length = 100)
  d1 <- data[,1]; d2 <- data[,2]
  sum1 <- sum(d1); slm1 <- sum(log(d1))
  sum2 <- sum(d2); slm2 <- sum(log(d2))
  nn <- nrow(data)
  nx <- length(x)
    
  mpmin <- function(a,b) {
    a[a > b] <- b[a > b]
    a
  }
  mpmax <- function(a,b) {
    a[a < b] <- b[a < b]
    a
  }
    
  if(method == "cfg") {
    if(!convex) {
      a <- numeric(nx)
      for(i in 1:nx)
        a[i] <- sum(log(mpmax((1-x[i]) * d1, x[i] * d2)))
      a <- (a - (1-x) * slm1 - x * slm2)/nn
      a <- pmin(1, pmax(exp(a), x, 1-x))
    }
    else {
      x2 <- seq(0, 1, length = 250)
      a <- numeric(250)
      for(i in 1:250)
        a[i] <- sum(log(mpmax((1-x2[i]) * d1, x2[i] * d2)))
      a <- (a - (1-x2) * slm1 - x2 * slm2)/nn
      a <- pmin(1, pmax(exp(a), x2, 1-x2))
      inch <- chull(x2, a)
      a <- a[inch] ; x2 <- x2[inch]
      a <- approx(x2, a, xout = x, method="linear")$y
    }
  }
   
  if(method == "pickands") {
    if(!convex) {
      a <- numeric(nx)
      if(madj == 2) {
        d1 <- d1/mean(d1)
        d2 <- d2/mean(d2)
      } 
      for(i in 1:nx)
        a[i] <- sum(mpmin(d1/x[i], d2/(1-x[i])))
      if(madj == 1)
        a <- a - x * sum1 - (1-x) * sum2 + nn
      a <- nn / a
      a <- pmin(1, pmax(a, x, 1-x))
    }
    else {
      x2 <- seq(0, 1, length = 250)
      a <- numeric(250)
      if(madj == 2) {
        d1 <- d1/mean(d1)
        d2 <- d2/mean(d2)
      } 
      for(i in 1:250)
        a[i] <- sum(mpmin(d1/x2[i], d2/(1-x2[i])))
      if(madj == 1)
        a <- a - x2 * sum1 - (1-x2) * sum2 + nn
      a <- nn / a
      a <- pmin(1, pmax(a, x2, 1-x2))
      inch <- chull(x2, a)
      a <- a[inch] ; x2 <- x2[inch]
      a <- approx(x2, a, xout = x, method="linear")$y
    }
  }   

  # Undocumented method: Tiago de Oliveira (1997)
  if(method == "tdo") {
    if(!convex) {
      a <- numeric(nx)
      for(i in 1:nx)
        a[i] <- sum(mpmin(x[i]/(1 + nn*d1), (1-x[i])/(1 + nn*d2)))
      a <- 1 - a/(1 + log(nn))
      a <- pmin(1, pmax(a, x, 1-x))
    }
    else {
      x2 <- seq(0, 1, length = 250)
      a <- numeric(250)
      for(i in 1:250)
        a[i] <- sum(mpmin(x2[i]/(1 + nn*d1), (1-x2[i])/(1 + nn*d2)))
      a <- 1 - a/(1 + log(nn))
      a <- pmin(1, pmax(a, x2, 1-x2))
      inch <- chull(x2, a)
      a <- a[inch] ; x2 <- x2[inch]
      a <- approx(x2, a, xout = x, method="linear")$y
    }
  }

  # Undocumented pot method: Beirlant et al (2004)
  if(method == "pot") {
    a <- numeric(nx)
    rr <- rowSums(1/data)
    rrk <- sort(rr, decreasing = TRUE)[k+1]
    for(i in 1:nx)
      a[i] <- sum(mpmax(x[i]/(d1 * rr), (1 - x[i])/(d2 * rr))[rr > rrk])
    a <- 2/k * a
    a0 <- 2/k * sum(1/(d2 * rr)[rr > rrk])
    a1 <- 2/k * sum(1/(d1 * rr)[rr > rrk])
    a <- a + 1 - (1-x) * a0 - x * a1
    a <- pmin(1, pmax(a, x, 1-x))
  } 
 
  if(plot || add) {
    if(!add)  { 
      plot(x, a, type="n", xlab = xlab, ylab = ylab, xlim = xlim,
        ylim = ylim, ...) 
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)  
    }
    lines(x, a, lty = lty, lwd = lwd, col = col) 
    return(invisible(list(x = x, y = a)))
  } 
  a
}

"qcbvnonpar"<- 
function(p = seq(0.75, 0.95, 0.05), data, epmar = FALSE, nsloc1 = NULL,
         nsloc2 = NULL, mint = 1, method = c("cfg","pickands","tdo"),
         convex = FALSE, madj = 0, kmar = NULL, plot = FALSE, add = FALSE,
         lty = 1, lwd = 1, col = 1, xlim = range(data[,1], na.rm = TRUE),
         ylim = range(data[,2], na.rm = TRUE), xlab = colnames(data)[1],
         ylab = colnames(data)[2], ...)
{
    if(mode(p) != "numeric" || p <= 0 || p >= 1)
      stop("`p' must be a vector of probabilities")
    method <- match.arg(method)
    nxv <- 100
    x <- seq(0, 1, length = nxv)
    ax <- abvnonpar(x = x, data = data, epmar = epmar, nsloc1 = nsloc1,
      nsloc2 = nsloc2, method = method, convex = convex, madj = madj,
      kmar = kmar, plot = FALSE)
    np <- length(p)
    qct <- list()
    p <- p^mint
    if(add) {
      xlim <- par("usr")[1:2]
      ylim <- par("usr")[1:2]
      if(par("xlog")) xlim <- 10^xlim
      if(par("ylog")) ylim <- 10^ylim
    }
    for(i in 1:np) {
      qct[[i]] <- -cbind(x/ax * log(p[i]), (1-x)/ax * log(p[i]))
      if(epmar) {
        qct[[i]] <- cbind(quantile(data[,1], probs = exp(-qct[[i]][,1]),
          na.rm = TRUE), quantile(data[,2], probs = exp(-qct[[i]][,2]),
          na.rm = TRUE))
      }
      else {
        if(is.null(kmar)) {
          # Transform from exponential margins
          mle.m1 <- fitted(fgev(data[,1], nsloc = nsloc1, std.err = FALSE))
          mle.m2 <- fitted(fgev(data[,2], nsloc = nsloc2, std.err = FALSE))
          mle.m1 <- mle.m1[c("loc","scale","shape")]
          mle.m2 <- mle.m2[c("loc","scale","shape")]
          qct[[i]] <- mtransform(qct[[i]], list(mle.m1, mle.m2), inv = TRUE)
        }
        else {
          if(!is.null(nsloc1) || !is.null(nsloc2))
            warning("ignoring `nsloc1' and `nsloc2' arguments")
          qct[[i]] <- mtransform(qct[[i]], kmar, inv = TRUE)
        }
      }
      qct[[i]][1,1] <- 1.5 * xlim[2]
      qct[[i]][nxv,2] <- 1.5 * ylim[2]
    }
    if((!is.null(nsloc1) || !is.null(nsloc2)) && !epmar && is.null(kmar)) {
        data <- fbvevd(data, model = "log", dep = 1, nsloc1 = nsloc1,
          nsloc2 = nsloc2, std.err = FALSE)$tdata
    }
    if(plot) {
      plot(data, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
      for(i in 1:np) lines(qct[[i]], lty = lty, lwd = lwd, col = col)
      return(invisible(qct))
    }
    if(add) {
      for(i in 1:np) lines(qct[[i]], lty = lty, lwd = lwd, col = col)
      return(invisible(qct))
    }
    qct
}

"amvnonpar"<- 
function(x = rep(1/d,d), data, d = 3, epmar = FALSE, nsloc = NULL,
    madj = 0, kmar = NULL, plot = FALSE, col = heat.colors(12),
    blty = 0, grid = if(blty) 150 else 50, lower = 1/3, ord = 1:3,
    lab = as.character(1:3), lcex = 1)
{ 
  if(!plot) {
    if(is.vector(x)) x <- as.matrix(t(x))
    if(!is.matrix(x) || ncol(x) != d)
      stop("`x' must be a vector/matrix with `d' elements/columns")
    if(any(x < 0, na.rm = TRUE))
      stop("`x' must be non-negative")
    rs <- rowSums(x)
    if(any(rs <= 0, na.rm = TRUE))
      stop("row(s) of `x' must have a positive sum")
    if(max(abs(rs[!is.na(rs)] - 1)) > 1e-6)
      warning("row(s) of `x' will be rescaled")
    x <- x/rs
  }
  if(missing(data) || ncol(data) != d)
    stop("data must have `d' columns")
  if(plot) {
    if(d == 2) stop("use abvnonpar for bivariate plots")
    if(d >= 4) stop("cannot plot in high dimensions")
  }
  if(epmar) {
    # Empirical transform to exponential
    data <- apply(data, 2, rank, na.last = "keep")
    nasm <- apply(data, 2, function(x) sum(!is.na(x)))
    data <- data / rep(nasm+1, each = nrow(data))
    data <- -log(data)
    # End empirical transform
  }
  if(!epmar) {
    # Parametric gev transform to exponential
    if(!is.null(kmar)) {
      data <- mtransform(data, kmar)
    }
    else {
      if(!is.null(nsloc)) {
        nslocmat <- list()
        if(!is.list(nsloc)) nsloc <- rep(list(nsloc), d)
        if(length(nsloc) != d)
          stop("`nsloc' should have `d' elements")
        for(k in 1:d) {
          if(!is.null(nsloc[[k]])) {
            if(is.vector(nsloc[[k]]))
              nsloc[[k]] <- data.frame(trend = nsloc[[k]])
            if(nrow(nsloc[[k]]) != nrow(data))
              stop("`nsloc' and data are not compatible")
            nslocmat[[k]] <- cbind(1, as.matrix(nsloc[[k]]))
          }
        }
      }   
      mles <- list()
      for(k in 1:d) {
        mle.m <- fitted(fgev(data[,k], nsloc = nsloc[[k]], std.err = FALSE))
        loc.mle.m <- mle.m[grep("^loc", names(mle.m))]
        if(is.null(nsloc[[k]])) loc.mle.m <- rep(loc.mle.m, nrow(data))
        else loc.mle.m <- nslocmat[[k]] %*% loc.mle.m
        mles[[k]] <- cbind(loc.mle.m, mle.m["scale"], mle.m["shape"])
      }
      data <- mtransform(data, mles)
    }
    # End parametric gev transform
  }
  
  data <- na.omit(data)
  depfn <- function(x, data, madj)
  {
    # quicker apply(am, 2, min)
    mpmin <- function(am,nr) {
      a <- am[1,]
      for(i in 2:nr)
        a[a > am[i,]] <- am[i, a > am[i,]]
      a
    }
    nn <- nrow(data)
    nx <- nrow(x)
    csum <- colSums(data)
    a <- numeric(nx)
    if(madj == 2)
      data <- nn * sweep(data, 2, csum, "/")
    for(i in 1:nx) {
      a[i] <- sum(mpmin(t(data)/x[i,], d))
    }
    if(madj == 1)
      a <- a - colSums(t(x) * csum) + nn
    a <- nn / a
    xrmax <- apply(x, 1, max)
    pmin(1, pmax(a, xrmax))
  }

  if(plot) {
    mz <- tvdepfn(depfn = depfn, col = col, blty = blty, grid = grid,
      lower = lower, ord = ord, lab = lab, lcex = lcex, data = data,
      madj = madj)
    return(invisible(mz))
  }
  depfn(x = x, data = data, madj = madj)
}



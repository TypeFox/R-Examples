
"evmc" <-
function(n, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
    margins = c("uniform","rweibull","frechet","gumbel"))
{
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct", "amix")
  m2 <- c(m1, c("log", "hr", "neglog"))
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy))
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta))
    warning("ignoring `beta' argument")

  nn <- as.integer(1)
  if(!(model %in% m1)) {
    if(length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
        stop("invalid argument for `dep'")
    if((model %in% c("log", "alog")) && dep > 1)
        stop("`dep' must be in the interval (0,1]")
    dep <- as.double(dep)
  }
  if(!(model %in% m2)) {
    if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
       max(asy) > 1) stop("invalid argument for `asy'")
    if(model == "alog" && (dep == 1 || any(asy == 0))) {
        asy <- c(0,0)
        dep <- 1
    }
    asy <- as.double(asy[c(2,1)])
  }
  if(!(model %in% m3)) {
    if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for `alpha'")
    if(length(beta) != 1 || mode(beta) != "numeric")
        stop("invalid argument for `beta'")
    if(model != "amix" && (alpha <= 0 || beta <= 0))
        stop("`alpha' and `beta' must be positive")
    if(model == "bilog" && any(c(alpha,beta) >= 1))
        stop("`alpha' and `beta' must be in the open interval (0,1)")
    if(model == "amix") {
      if(alpha < 0)
        stop("`alpha' must be non-negative")
      if((alpha + beta) > 1)
        stop("`alpha' + `beta' cannot be greater than one")
      if((alpha + 2*beta) > 1)
        stop("`alpha' + `2*beta' cannot be greater than one")
      if((alpha + 3*beta) < 0)
        stop("`alpha' + `3*beta' must be non-negative")
      alpha <- as.double(alpha + 3*beta)
      beta <- as.double(-beta)
    }
    else {
      alpha <- as.double(beta)
      beta <- as.double(alpha)
    }
  }
    
  evmc <- runif(n)
  for(i in 2:n) {
    evmc[c(i,i-1)] <- switch(model,
      log = .C("rbvlog", nn, dep, sim = evmc[c(i,i-1)], PACKAGE = "evd")$sim,
      alog = .C("rbvalog", nn, dep, asy, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      hr = .C("rbvhr", nn, dep, sim = evmc[c(i,i-1)], PACKAGE = "evd")$sim,
      neglog = .C("rbvneglog", nn, dep, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      aneglog = .C("rbvaneglog", nn, dep, asy, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      bilog = .C("rbvbilog", nn, alpha, beta, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      negbilog = .C("rbvnegbilog", nn, alpha, beta, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      ct = .C("rbvct", nn, alpha, beta, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim,
      amix = .C("rbvamix", nn, alpha, beta, sim = evmc[c(i,i-1)],
        PACKAGE = "evd")$sim)
  }

  switch(match.arg(margins),
    frechet = -1/log(evmc), uniform = evmc,
    rweibull = log(evmc), gumbel = -log(-log(evmc))) 
}

"marma" <-
function(n, p = 0, q = 0, psi, theta, init = rep(0, p), n.start = p,
    rand.gen = rfrechet, ...)
{
    if(missing(psi)) psi <- numeric(0)
    if(missing(theta)) theta <- numeric(0)
    if(length(psi) != p || !is.numeric(psi) || any(psi < 0))
      stop("`par' must be a non-negative vector of length `p'")
    if(length(theta) != q || !is.numeric(theta) || any(theta < 0))
      stop("`theta' must be a non-negative vector of length `q'")
    if(length(init) != p || !is.numeric(init) || any(init < 0))
      stop("`init' must be a non-negative vector of length `p'")
    
    marma <- c(init, numeric(n + n.start - p))
    theta <- c(1, theta)
    innov <- rand.gen(n.start + n + q, ...)
    for(i in 1:(n + n.start - p))
      marma[i+p] <- max(c(psi * marma[(i+p-1):i], theta * innov[(i+q):i]))
    if(n.start) marma <- marma[-(1:n.start)]
    marma
}

"mma" <-
function(n, q = 1, theta, rand.gen = rfrechet, ...)
{
    marma(n = n, q = q, theta = theta, rand.gen = rand.gen, ...)
}

"mar" <-
function(n, p = 1, psi, init = rep(0, p), n.start = p, rand.gen =
         rfrechet, ...)
{
    marma(n = n, p = p, psi = psi, init = init, n.start = n.start,
          rand.gen = rand.gen, ...)
}

"clusters"<-
function(data, u, r = 1, ulow = -Inf, rlow = 1, cmax = FALSE, keep.names = TRUE,
    plot = FALSE, xdata = seq(along = data), lvals = TRUE, lty = 1, lwd = 1,
    pch = par("pch"), col = if(n > 250) NULL else "grey", xlab = "Index",
    ylab = "Data", ...)
{
    n <- length(data)
    if(length(u) != 1) u <- rep(u, length.out = n)
    if(length(ulow) != 1) ulow <- rep(ulow, length.out = n)
    if(any(ulow > u)) stop("`u' cannot be less than `ulow'")
    if(is.null(names(data)) && keep.names) names(data) <- 1:n
    if(!keep.names) names(data) <- NULL
    high <- as.double((data > u) & !is.na(data))
    high2 <- as.double((data > ulow) | is.na(data))
    clstrs <- .C("clusters", high, high2, n, as.integer(r),
        as.integer(rlow), clstrs = double(3*n), PACKAGE = "evd")$clstrs
    clstrs <- matrix(clstrs, n, 3)
    start <- clstrs[,2] ; end <- clstrs[,3]
    splvec <- clstrs[,1]
    start <- as.logical(start)
    end <- as.logical(end)
    clstrs <- split(data, splvec)
    names(clstrs) <- paste("cluster", names(clstrs), sep = "")
    if(any(!splvec)) clstrs <- clstrs[-1]
    nclust <- length(clstrs)
    acs <- sum(high)/nclust
    if(plot) {
      if(length(xdata) != length(data))
        stop("`xdata' and `data' have different lengths")
      if(any(is.na(xdata)))
        stop("`xdata' cannot contain missing values")
      if(any(duplicated(xdata)))
        stop("`xdata' cannot contain duplicated values")
      eps <- min(diff(xdata))/2
      start <- xdata[start] - eps
      end <- xdata[end] + eps
      plot(xdata, data, xlab = xlab, ylab = ylab, type = "n", ...)
      if(!is.null(col) && nclust > 0.5) {
        for(i in 1:nclust) {
          xvl <- c(start[i], end[i], end[i], start[i])
          polygon(xvl, rep(par("usr")[3:4], each = 2), col = col)
        }
      }
      if(length(u) == 1) abline(h = u, lty = lty, lwd = lwd)
      else lines(xdata, u, lty = lty, lwd = lwd)
      if(lvals) {
        if(length(ulow) == 1) abline(h = ulow, lty = lty, lwd = lwd)
        else lines(xdata, ulow, lty = lty, lwd = lwd)
      }
      else {
        high <- as.logical(high)
        xdata <- xdata[high]
        data <- data[high]
      }
      points(xdata, data, pch = pch)
    }
    if(cmax) {
      if(keep.names)
        nmcl <- unlist(lapply(clstrs, function(x) names(x)[which.max(x)]))
      clstrs <- as.numeric(unlist(lapply(clstrs, max, na.rm = TRUE)))
      if(keep.names) names(clstrs) <- nmcl 
    }
    attributes(clstrs)$acs <- acs
    if(plot) return(invisible(clstrs))   
    clstrs
}

"exi"<-
function (data, u, r = 1, ulow = -Inf, rlow = 1) 
{
    n <- length(data)
    if (length(u) != 1) 
        u <- rep(u, length.out = n)
    if (length(ulow) != 1) 
        ulow <- rep(ulow, length.out = n)
    if (any(ulow > u)) 
        stop("`u' cannot be less than `ulow'")

    if(r > 0.5) {
      clstrs <- clusters(data, u = u, r = r, ulow = ulow, 
        rlow = rlow, keep.names = FALSE)
      exindex <- 1/attributes(clstrs)$acs
    }
    else {
      extms <- which(data > u)
      nn <- length(extms)
      if(nn == 0) return(NaN)
      if(nn == 1) return(1)
      iextms <- extms[-1] - extms[-nn]
      if(max(iextms) > 2.5) {
        den <- log(nn - 1) + log(sum((iextms - 1) * (iextms - 2)))
        exindex <- log(2) + 2*log(sum(iextms - 1)) - den
        exindex <- min(1, exp(exindex))
      }
      else {
        den <- log(nn - 1) + log(sum(iextms^2))
        exindex <- log(2) + 2*log(sum(iextms)) - den
        exindex <- min(1, exp(exindex))
      }
    }
    exindex
}

exiplot <-
function (data, tlim, r = 1, ulow = -Inf, rlow = 1, add = FALSE, 
    nt = 100, lty = 1, xlab = "Threshold", ylab = "Ext. Index",
    ylim = c(0,1), ...) 
{
    nn <- length(data)
    if (all(data <= tlim[2])) 
        stop("upper limit for threshold is too high")
    u <- seq(tlim[1], tlim[2], length = nt)
    x <- numeric(nt)
    for (i in 1:nt) {
        x[i] <- exi(data, u = u[i], r = r, ulow = ulow, rlow = rlow)
    }
    if(add) {
      lines(u, x, lty = lty, ...)
    } 
    else {
      plot(u, x, type = "l", lty = lty, xlab = xlab, 
        ylab = ylab, ylim = ylim, ...)
    }
    invisible(list(x = u, y = x))
}


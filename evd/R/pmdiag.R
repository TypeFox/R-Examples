
"mrlplot"<-
function(data, tlim, pscale = FALSE, nt = max(100, length(data)), lty =
    c(2,1,2), col = 1, conf = 0.95, main = "Mean Residual Life Plot",
    xlab = "Threshold", ylab = "Mean Excess", ...)
{
    data <- sort(data[!is.na(data)])
    nn <- length(data)
    if(nn <= 5) stop("`data' has too few non-missing values")
    if(missing(tlim)) {
      tlim <- c(data[1], data[nn - 4])
      tlim <- tlim - .Machine$double.eps^0.5
    }
    if(all(data <= tlim[2]))
      stop("upper limit for threshold is too high")
    u <- seq(tlim[1], tlim[2], length = nt)
    if(pscale) { 
      tlim[1] <- mean(data <= tlim[1], na.rm = TRUE)
      tlim[2] <- mean(data <= tlim[2], na.rm = TRUE)
      pvec <- seq(tlim[1], tlim[2], length = nt)
      u <- quantile(data, probs = pvec, na.rm = TRUE)
    }
    x <- matrix(NA, nrow = nt, ncol = 3, dimnames = list(NULL,
        c("lower", "mrl", "upper")))
    for(i in 1:nt) {
        data <- data[data > u[i]]
	x[i,2] <- mean(data - u[i])
	sdev <- sqrt(var(data))
        sdev <- (qnorm((1 + conf)/2) * sdev)/sqrt(length(data))
	x[i,1] <- x[i,2] - sdev
	x[i,3] <- x[i,2] + sdev
    }
    if(pscale) {
      u <- pvec
      if(missing(xlab)) xlab <- "Threshold probability"
    }
    matplot(u, x, type = "l", lty = lty, col = col, main = main,
            xlab = xlab, ylab = ylab, ...)
    invisible(list(x = u, y = x))
}

"tcplot"<-
function(data, tlim, model = c("gpd", "pp"), pscale = FALSE, cmax = FALSE,
    r = 1, ulow = -Inf, rlow = 1, nt = 25, which = 1:npar, conf = 0.95,
    lty = 1, lwd = 1, type = "b", cilty = 1, vci = TRUE, xlab, xlim, ylabs,
    ylims, ask = nb.fig < length(which) && dev.interactive(), ...)
{
    model <- match.arg(model)
    u <- seq(tlim[1], tlim[2], length = nt)
    if(pscale) { 
      tlim[1] <- mean(data <= tlim[1], na.rm = TRUE)
      tlim[2] <- mean(data <= tlim[2], na.rm = TRUE)
      pvec <- seq(tlim[1], tlim[2], length = nt)
      u <- quantile(data, probs = pvec, na.rm = TRUE)
    }
    locs <- scls <- shps <- matrix(NA, nrow = nt, ncol = 3)
    dimnames(locs) <- list(round(u,2), c("lower", "loc", "upper"))
    dimnames(shps) <- list(round(u,2), c("lower", "shape", "upper"))
    if(model == "gpd") {
      pname <- "mscale"
      npar <- 2
    }
    if(model == "pp") {
      pname <- "scale"
      npar <- 3
    }
    dimnames(scls) <- list(round(u,2), c("lower", pname, "upper"))
    z <- fpot(data, u[1], model = model, cmax = cmax, r = r, ulow = ulow,
              rlow = rlow, corr = TRUE, ...)
    stvals <- as.list(round(fitted(z), 3))
    for(i in 1:nt) {
      z <- fpot(data, u[i], model = model, start = stvals, cmax = cmax,
                r = r, ulow = ulow, rlow = rlow, corr = TRUE, ...)
      stvals <- as.list(fitted(z))
      mles <- fitted(z)
      stderrs <- std.errors(z)
      cnst <- qnorm((1 + conf)/2)
      shp <- mles["shape"]
      scl <- mles["scale"]
      shpse <- stderrs["shape"]
      sclse <- stderrs["scale"]
      if(model == "pp") {
        loc <- mles["loc"]  
        locse <- stderrs["loc"]
        locs[i,] <- c(loc - cnst*locse, loc, loc + cnst*locse)
      } 
      if(model == "gpd") {
        scl <- scl - shp*u[i]
        covar <- z$corr[1,2] * prod(stderrs)
        sclse <- sqrt(sclse^2 - 2*u[i]*covar + (u[i]*shpse)^2)
      }
      scls[i,] <- c(scl - cnst*sclse, scl, scl + cnst*sclse)
      shps[i,] <- c(shp - cnst*shpse, shp, shp + cnst*shpse)      
   }
   show <- rep(FALSE, npar)
   show[which] <- TRUE
   nb.fig <- prod(par("mfcol"))
   if (ask) {
     op <- par(ask = TRUE)
     on.exit(par(op))
   }
   if(pscale) u <- pvec
   if(missing(xlim)) xlim <- tlim
   if(missing(xlab)) {
     xlab <- "Threshold"
     if(pscale) xlab <- "Threshold probability"
   }
   if(model == "pp") {
     ylab <- c("Location","Scale","Shape")
     if(!missing(ylabs)) ylab[show] <- ylabs
     ylim <- rbind(range(locs), range(scls), range(shps))
     if(!missing(ylims)) ylim[show,] <- ylims
     if(show[1]) {
       matplot(u, locs, type = "n", xlab = xlab, ylab = ylab[1],
            xlim = xlim, ylim = ylim[1,])
       lines(u, locs[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, locs[,1], u, locs[,3], lty = cilty)
       else {
         lines(u, locs[,1], lty = cilty)
         lines(u, locs[,3], lty = cilty)
       }
     }
     if(show[2]) {
       matplot(u, scls, type = "n", xlab = xlab, ylab = ylab[2],
            xlim = xlim, ylim = ylim[2,])
       lines(u, scls[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, scls[,1], u, scls[,3], lty = cilty)
       else {
         lines(u, scls[,1], lty = cilty)
         lines(u, scls[,3], lty = cilty)
       }
     }
     if(show[3]) {
       matplot(u, shps, type = "n", xlab = xlab, ylab = ylab[3],
            xlim = xlim, ylim = ylim[3,])
       lines(u, shps[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, shps[,1], u, shps[,3], lty = cilty)
       else {
         lines(u, shps[,1], lty = cilty)
         lines(u, shps[,3], lty = cilty)
       }
     }
     rtlist <- list(locs = locs, scales = scls, shapes = shps) 
   }
   if(model == "gpd") {
     ylab <- c("Modified Scale","Shape")
     if(!missing(ylabs)) ylab[show] <- ylabs
     ylim <- rbind(range(scls), range(shps))
     if(!missing(ylims)) ylim[show,] <- ylims
     if(show[1]) {
       matplot(u, scls, type = "n", xlab = xlab, ylab = ylab[1],
            xlim = xlim, ylim = ylim[1,])
       lines(u, scls[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, scls[,1], u, scls[,3], lty = cilty)
       else {
         lines(u, scls[,1], lty = cilty)
         lines(u, scls[,3], lty = cilty)
       }
     }
     if(show[2]) {
       matplot(u, shps, type = "n", xlab = xlab, ylab = ylab[2],
            xlim = xlim, ylim = ylim[2,])
       lines(u, shps[,2], lty = lty, lwd = lwd, type = type)
       if(vci) segments(u, shps[,1], u, shps[,3], lty = cilty)
       else {
         lines(u, shps[,1], lty = cilty)
         lines(u, shps[,3], lty = cilty)
       }
     }
     rtlist <- list(scales = scls, shapes = shps) 
   } 
   invisible(rtlist)
}

"chiplot"<- 
function(data, nq = 100, qlim = NULL, which = 1:2, conf = 0.95, trunc = TRUE, spcases = FALSE, lty = 1, cilty = 2, col = 1, cicol = 1, xlim = c(0,1), ylim1 = c(-1,1), ylim2 = c(-1,1), main1 = "Chi Plot", main2 = "Chi Bar Plot", xlab = "Quantile", ylab1 = "Chi", ylab2 = "Chi Bar", ask = nb.fig < length(which) && dev.interactive(), ...)
{
    data <- na.omit(data)
    n <- nrow(data)
    data <- cbind(rank(data[, 1])/(n + 1), rank(data[, 2])/(n + 1))
    rowmax <- apply(data, 1, max)
    rowmin <- apply(data, 1, min)
    eps <- .Machine$double.eps^0.5
    qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)
    if(!is.null(qlim)) {
      if(qlim[1] < qlim2[1]) stop("lower quantile limit is too low")
      if(qlim[2] > qlim2[2]) stop("upper quantile limit is too high")
      if(qlim[1] > qlim[2]) stop("lower quantile limit is less than upper quantile limit")
    } else qlim <- qlim2
    u <- seq(qlim[1], qlim[2], length = nq)

    cu <- cbaru <- numeric(nq)
    for(i in 1:nq) cu[i] <- mean(rowmax < u[i])
    for(i in 1:nq) cbaru[i] <- mean(rowmin > u[i])
    chiu <- 2 - log(cu)/log(u)
    chibaru <- (2 * log(1 - u))/log(cbaru) - 1

    cnst <- qnorm((1 + conf)/2)
    varchi <- ((1/log(u)^2 * 1)/cu^2 * cu * (1 - cu))/n
    varchi <- cnst * sqrt(varchi)
    varchibar <- (((4 * log(1 - u)^2)/(log(cbaru)^4 * cbaru^2)) * cbaru * (
      1 - cbaru))/n
    varchibar <- cnst * sqrt(varchibar)
    chiu <- cbind(chilow = chiu-varchi, chi = chiu, chiupp = chiu+varchi) 
    chibaru <- cbind(chiblow = chibaru-varchibar, chib = chibaru, chibupp =
      chibaru+varchibar)

    chiulb <- 2-log(pmax(2*u-1,0))/log(u)
    chibarulb <- 2*log(1-u)/log(1-2*u+pmax(2*u-1,0)) - 1
    if(trunc) {
      chiu[chiu > 1] <- 1
      chibaru[chibaru > 1] <- 1
      chiu <- apply(chiu, 2, function(x) pmax(x, chiulb))
      chibaru <- apply(chibaru, 2, function(x) pmax(x, chibarulb))
    }

    show <- logical(2)
    show[which] <- TRUE
    lty <- c(cilty, lty, cilty)
    col <- c(cicol, col, cicol)
    nb.fig <- prod(par("mfcol"))
    if (ask) {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }

    if(show[1]) {
      matplot(u, chiu, type = "l", lty = lty, col = col, xlim = xlim, ylim =
         ylim1, main = main1, xlab = xlab, ylab = ylab1, ...)
      if(spcases) {
        segments(qlim[1],0,qlim[2],0, lty = 5, col = "grey")
        segments(qlim[1],1,qlim[2],1, lty = 5, col = "grey")
        lines(u, chiulb, lty = 5, col = "grey")
      }
    }

    if(show[2]) {
      matplot(u, chibaru, type = "l", lty = lty, col = col, xlim = xlim,
         ylim = ylim2, main = main2, xlab = xlab, ylab = ylab2, ...)
      if(spcases) {
        segments(qlim[1],0,qlim[2],0, lty = 5, col = "grey")
        segments(qlim[1],1,qlim[2],1, lty = 5, col = "grey")
        lines(u, chibarulb, lty = 5, col = "grey")
      }
    }

    plvals <- list(quantile = u, chi = chiu, chibar = chibaru)
    if(!show[1]) plvals$chi <- NULL
    if(!show[2]) plvals$chib <- NULL
    invisible(plvals)    
}

## Bivariate Threshold Choice ##

"bvtcplot"<- 
function(x, spectral = FALSE, xlab, ylab, ...)
{
  if(!is.matrix(x) && !is.data.frame(x))
    stop("`x' must be a matrix or data frame")
  if(ncol(x) != 2)
	  stop("`x' has incorrect number of columns")
	
	x <- x[complete.cases(x),]
	nn <- nrow(x)
	ula <- apply(x, 2, rank)/(nn + 1)
  fla <- -1/log(ula)
  rr <- rowSums(fla); ww <- fla/rr
	
  rro <- sort(rr, decreasing = TRUE)[-1]
  k <- 1:(nn-1)
	k0 <- max(which(rro*k/nn > 2))
  if(!spectral) {
	  if(missing(xlab)) xlab <- "k"
		if(missing(ylab)) ylab <- "H([0,1])"
	  plot(k, rro*k/nn, xlab = xlab, ylab = ylab, ...)
	  abline(h = 2, v = k0)
		return(invisible(list(x = k, y = rro*k/nn, k0 = k0)))
	}
	
  xx <- yy <- seq(0, 1, len = 100)
  for(k in 1:100) yy[k] <- sum(rr > rro[k0] & ww[,1] <= xx[k])
	if(missing(xlab)) xlab <- "w"
	if(missing(ylab)) ylab <- "H([0,w])"
  plot(xx, 2/k0 * yy, type = "l", xlab = xlab, ylab = ylab, ...)
  abline(h = c(0,2))
	return(invisible(list(x = xx, y = 2/k0 * yy, k0 = k0)))
}

## Hypothesis test for independence ##

"evind.test"<-
function(x, method = c("ratio", "score"), verbose = FALSE)
{
  method <- match.arg(method)
  if(!is.matrix(x) && !is.data.frame(x))
    stop("`x' must be a matrix or data frame")
  if(ncol(x) != 2)
	  stop("`x' has incorrect number of columns")
	dname <- paste(deparse(substitute(x)))
		
	if(method == "ratio") {
		meth <- "Likelihood Ratio Test Of Independence"
	  fobj1 <- fbvevd(x, model = "log")
		estimate <- fitted(fobj1)
		if(!verbose) estimate <- estimate["dep"]
	  fobj2 <- fbvevd(x, model = "log", dep = 1)
	  lrt <- anova(fobj1, fobj2, half = TRUE)
		stat <- c(norm.llhratio = lrt[["Chisq"]][2])
		pval <- c(p.value = lrt[["Pr(>chisq)"]][2])
 }
	if(method == "score") {
		meth <- "Score Test Of Independence"
		fobj1 <- fbvevd(x, model = "log")
		estimate <- fitted(fobj1)
		if(!verbose) estimate <- estimate["dep"]
		fobj2 <- fbvevd(x, model = "log", dep = 1)
		ft <- fitted(fobj2)
	  mmles <- list(ft[c("loc1","scale1","shape1")], 
			ft[c("loc2","scale2","shape2")])
    xtr <- mtransform(x, mmles)
		xtr <- xtr[complete.cases(xtr),]
		nn <- nrow(xtr)
    rsm <- rowSums(xtr)
    rsl <- rowSums(xtr * log(xtr))
    tawn <- rsl - log(apply(xtr, 1, prod)) - (rsm - 2) * log(rsm) - 1/rsm
		tawn <- (nn/2 * log(nn))^(-1/2) * sum(tawn)
		stat <- c(norm.score = tawn)
		pval <- c(p.value = pnorm(tawn))
	}
	rval <- list(statistic = stat, p.value = pval, estimate = estimate, null.value = 
	  c(dependence = "independence"), alternative = "greater", method = meth, 
		data.name = dname)
  class(rval) <- "htest"
  return(rval)
}




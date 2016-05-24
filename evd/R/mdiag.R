### Univariate GEV and POT Models ###

"plot.uvevd" <-  function(x, which = 1:4, main, ask = nb.fig <
     length(which) && dev.interactive(), ci = TRUE, cilwd = 1,
     adjust = 1, jitter = FALSE, nplty = 2, ...) 
{
    if (!inherits(x, "uvevd")) 
        stop("Use only with `'uvevd objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")
    nb.fig <- prod(par("mfcol"))
    show <- rep(FALSE, 4)
    show[which] <- TRUE
   if(missing(main)) {
      main <- c("Probability Plot", "Quantile Plot", "Density Plot",
        "Return Level Plot")
    }
    else {
      if(length(main) != length(which))
        stop("number of plot titles is not correct")
      main2 <- character(4)
      main2[show] <- main
      main <- main2
    }
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        pp(x, ci = ci, cilwd = cilwd, main = main[1], xlim = c(0,1),
           ylim = c(0,1), ...)
    }
    if (show[2]) {
        qq(x, ci = ci, cilwd = cilwd, main = main[2], ...)
    }
    if (show[3]) {
        dens(x, adjust = adjust, nplty = nplty, jitter = jitter,
             main = main[3], ...)
    }
    if (show[4]) {
        rl(x, ci = ci, cilwd = cilwd, main = main[4], ...)
    }
    invisible(x)
}

"qq" <- function (x, ...) UseMethod("qq")
"pp" <- function (x, ...) UseMethod("pp")
"rl" <- function (x, ...) UseMethod("rl")
"dens" <- function (x, ...) UseMethod("dens")

"qq.gev" <-  function(x, ci = TRUE, cilwd = 1, main = "Quantile Plot", xlab = "Model", ylab = "Empirical", ...)
{
    quant <- qgev(ppoints(x$tdata), loc = x$loc,
                 scale = x$param["scale"], shape = x$param["shape"])
    if(!ci) {
      plot(quant, sort(x$tdata), main = main, xlab = xlab, ylab = ylab, ...)
      abline(0, 1)
    }
    else {
      samp <- rgev(x$n*99, loc = x$loc,
                 scale = x$param["scale"], shape = x$param["shape"])
      samp <- matrix(samp, x$n, 99)
      samp <- apply(samp, 2, sort)
      samp <- apply(samp, 1, sort)
      env <- t(samp[c(3,97),])
      rs <- sort(x$tdata)
      matplot(quant, cbind(rs,env), main = main, xlab = xlab, ylab = ylab,
              type = "pnn", pch = 4, ...)
      xyuser <- par("usr")
      smidge <- min(diff(c(xyuser[1], quant, xyuser[2])))/2
      smidge <- max(smidge, (xyuser[2] - xyuser[1])/200)
      segments(quant-smidge, env[,1], quant+smidge, env[,1], lwd = cilwd)
      segments(quant-smidge, env[,2], quant+smidge, env[,2], lwd = cilwd)
      abline(0, 1)
    }
    invisible(list(x = quant, y = sort(x$tdata)))
}

"pp.gev" <-  function(x, ci = TRUE, cilwd = 1, main = "Probability Plot", xlab = "Empirical", ylab = "Model", ...)
{
    ppx <- ppoints(x$n)
    probs <- pgev(sort(x$tdata), loc = x$loc,
                 scale = x$param["scale"], shape = x$param["shape"])
    if(!ci) {
        plot(ppx, probs, main = main, xlab = xlab, ylab = ylab, ...)
        abline(0, 1)
    }
    else {
        samp <- rgev(x$n*99, loc = x$loc,
                   scale = x$param["scale"], shape = x$param["shape"])
        samp <- matrix(samp, x$n, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        env[,1] <- pgev(env[,1], loc = x$loc,
                    scale = x$param["scale"], shape = x$param["shape"])
        env[,2] <- pgev(env[,2], loc = x$loc,
                    scale = x$param["scale"], shape = x$param["shape"])
        matplot(ppx, cbind(probs, env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, ...)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], ppx, xyuser[2])))/2
        smidge <- max(smidge, (xyuser[2] - xyuser[1])/200)
        segments(ppx-smidge, env[,1], ppx+smidge, env[,1], lwd = cilwd)
        segments(ppx-smidge, env[,2], ppx+smidge, env[,2], lwd = cilwd)
        abline(0, 1)
    }
    invisible(list(x = ppoints(x$n), y = probs))
}

"rl.gev" <-  function(x, ci = TRUE, cilwd = 1, main = "Return Level Plot", xlab = "Return Period", ylab = "Return Level", ...)
{
    ppx <- ppoints(x$tdata)
    rps <- c(1.001,10^(seq(0,3,len=200))[-1])
    p.upper <- 1/rps
    rlev <- qgev(p.upper, loc = x$loc, scale = x$param["scale"],
              shape = x$param["shape"], lower.tail = FALSE)
    if(!ci) {
        plot(-1/log(ppx), sort(x$tdata),log = "x", main = main,
             xlab = xlab, ylab = ylab, ...)
        lines(-1/log(1-p.upper), rlev)
    }
    else {
        samp <- rgev(x$n*99, loc = x$loc,
                   scale = x$param["scale"], shape = x$param["shape"])
        samp <- matrix(samp, x$n, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        rs <- sort(x$tdata)
        matplot(-1/log(ppx), cbind(rs,env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, log = "x", ...)
        lines(-1/log(1-p.upper), rlev)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], log10(-1/log(ppx)), xyuser[2])))/2
        smidge <- max(smidge, (xyuser[2] - xyuser[1])/200)
        segments((-1/log(ppx))*exp(-smidge), env[,1],
                 (-1/log(ppx))*exp(smidge), env[,1], lwd = cilwd)
        segments((-1/log(ppx))*exp(-smidge), env[,2],
                 (-1/log(ppx))*exp(smidge), env[,2], lwd = cilwd)
    }
    invisible(list(x = -1/log(1-p.upper), y = rlev))
}

"dens.gev" <-  function(x, adjust = 1, nplty = 2, jitter = FALSE, main = "Density Plot", xlab = "Quantile", ylab = "Density", ...)
{
    xlimit <- range(x$tdata)
    xlimit[1] <- xlimit[1] - diff(xlimit) * 0.075
    xlimit[2] <- xlimit[2] + diff(xlimit) * 0.075
    xvec <- seq(xlimit[1], xlimit[2], length = 100)
    dens <- dgev(xvec, loc = x$loc, scale = x$param["scale"],
                shape = x$param["shape"])
    plot(spline(xvec, dens), main = main, xlab = xlab, ylab = ylab,
         type = "l", ...)
    if(jitter) rug(jitter(x$tdata))
    else rug(x$tdata)
    lines(density(x$tdata, adjust = adjust), lty = nplty)
    invisible(list(x = xvec, y = dens))
}

"qq.pot" <-  function(x, ci = TRUE, cilwd = 1, main = "Quantile Plot", xlab = "Model", ylab = "Empirical", ...)
{
    quant <- qgpd(ppoints(x$nhigh), loc = x$threshold,
                 scale = x$scale, shape = x$param["shape"])
    if(!ci) {
      plot(quant, sort(x$exceedances), main = main, xlab = xlab,
           ylab = ylab, ...)
      abline(0, 1)
    }
    else {
      samp <- rgpd(x$nhigh*99, loc = x$threshold,
                 scale = x$scale, shape = x$param["shape"])
      samp <- matrix(samp, x$nhigh, 99)
      samp <- apply(samp, 2, sort)
      samp <- apply(samp, 1, sort)
      env <- t(samp[c(3,97),])
      rs <- sort(x$exceedances)
      matplot(quant, cbind(rs,env), main = main, xlab = xlab, ylab = ylab,
              type = "pnn", pch = 4, ...)
      xyuser <- par("usr")
      smidge <- min(diff(c(xyuser[1], quant, xyuser[2])))/2
      smidge <- max(smidge, (xyuser[2] - xyuser[1])/200)
      segments(quant-smidge, env[,1], quant+smidge, env[,1], lwd = cilwd)
      segments(quant-smidge, env[,2], quant+smidge, env[,2], lwd = cilwd)
      abline(0, 1)
    }
    invisible(list(x = quant, y = sort(x$exceedances)))
}

"pp.pot" <-  function(x, ci = TRUE, cilwd = 1, main = "Probability Plot", xlab = "Empirical", ylab = "Model", ...)
{
    ppx <- ppoints(x$nhigh)
    probs <- pgpd(sort(x$exceedances), loc = x$threshold,
                 scale = x$scale, shape = x$param["shape"])
    if(!ci) {
        plot(ppx, probs, main = main, xlab = xlab, ylab = ylab, ...)
        abline(0, 1)
    }
    else {
        samp <- rgpd(x$nhigh*99, loc = x$threshold,
                   scale = x$scale, shape = x$param["shape"])
        samp <- matrix(samp, x$nhigh, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        env[,1] <- pgpd(env[,1], loc = x$threshold,
                    scale = x$scale, shape = x$param["shape"])
        env[,2] <- pgpd(env[,2], loc = x$threshold,
                    scale = x$scale, shape = x$param["shape"])
        matplot(ppx, cbind(probs, env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, ...)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], ppx, xyuser[2])))/2
        smidge <- max(smidge, (xyuser[2] - xyuser[1])/200)
        segments(ppx-smidge, env[,1], ppx+smidge, env[,1], lwd = cilwd)
        segments(ppx-smidge, env[,2], ppx+smidge, env[,2], lwd = cilwd)
        abline(0, 1)
    }
    invisible(list(x = ppoints(x$nhigh), y = probs))
}

"rl.pot" <- function(x, ci = TRUE, cilwd = 1, main = "Return Level Plot", xlab = "Return Period", ylab = "Return Level", ...)
{
    rpstmfc <- c(1.001,10^(seq(0,3,len=200))[-1])
    rlev <- qgpd(1/rpstmfc, loc = x$threshold, scale = x$scale,
              shape = x$param["shape"], lower.tail = FALSE)
    mfc <- x$npp * x$nhigh/length(x$data)
    rps <- rpstmfc/mfc
    ppx <- 1/(mfc * (1 - ppoints(x$nhigh)))
    if(!ci) {
        plot(ppx, sort(x$exceedances), log = "x", main =
             main, xlab = xlab, ylab = ylab, ...)
        lines(rps, rlev)
    }
    else {
        samp <- rgpd(x$nhigh*99, loc = x$threshold,
                   scale = x$scale, shape = x$param["shape"])
        samp <- matrix(samp, x$nhigh, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        rs <- sort(x$exceedances)
        matplot(ppx, cbind(rs,env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, log = "x", ...)
        lines(rps, rlev)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], log10(ppx), xyuser[2])))/2
        smidge <- max(smidge, (xyuser[2] - xyuser[1])/200)
        segments(ppx*exp(-smidge), env[,1], ppx*exp(smidge),
                 env[,1], lwd = cilwd)
        segments(ppx*exp(-smidge), env[,2], ppx*exp(smidge),
                 env[,2], lwd = cilwd)
    }
    invisible(list(x = rps, y = rlev))
}

"dens.pot" <-  function(x, adjust = 1, nplty = 2, jitter = FALSE, main = "Density Plot", xlab = "Quantile", ylab = "Density", ...)
{
    xlimit <- range(x$exceedances)
    xlimit[2] <- xlimit[2] + diff(xlimit) * 0.075
    xvec <- seq(xlimit[1], xlimit[2], length = 100)
    dens <- dgpd(xvec, loc = x$threshold, scale = x$scale,
                shape = x$param["shape"])
    plot(spline(xvec, dens), main = main, xlab = xlab, ylab = ylab,
         type = "l", ...)
    if(jitter) rug(jitter(x$exceedances))
    else rug(x$exceedances)
    flipexceed <- c(x$exceedances, 2*x$threshold - x$exceedances)
    flipdens <- density(flipexceed, adjust = adjust, from = xlimit[1],
        to = xlimit[2])
    flipdens$y <- 2*flipdens$y
    lines(flipdens, lty = nplty)
    invisible(list(x = xvec, y = dens))
}

### Bivariate EVD Models ###

"plot.bvevd" <-  function(x, mar = 0, which = 1:6, main,
     ask = nb.fig < length(which) && dev.interactive(), ci = TRUE,
     cilwd = 1, grid = 50, legend = TRUE, nplty = 2,
     blty = 3, method = "cfg", convex = FALSE, rev = FALSE,
     p = seq(0.75, 0.95, 0.05), mint = 1, half = FALSE, ...) 
{
    if (!inherits(x, "bvevd")) 
        stop("Use only with `bvevd' objects")
    nb.fig <- prod(par("mfcol"))
    if(mar == 1 || mar == 2) {
        indx <- paste(c("loc","scale","shape"), as.character(mar), sep="")
        tdata <- na.omit(x$tdata[, mar])
        n <- length(tdata)
        param <- x$param[indx]
        names(param) <- c("loc","scale","shape")
        gev.mar <- structure(list(param = param, tdata = tdata, n = n,
           loc = param["loc"]), class = c("gev", "uvevd", "evd"))
        if(missing(which)) which <- 1:4
        plot(gev.mar, which = which, main = main, ask = ask, ci = ci,
             cilwd = cilwd, ...)
        return(invisible(x))
    }
    if (!is.numeric(which) || any(which < 1) || any(which > 6)) 
        stop("`which' must be in 1:6")
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    if(missing(main)) {
      main <- c("Conditional Plot One", "Conditional Plot Two",
        "Density Plot", "Dependence Function", "Quantile Curves",
        "Spectral Density")
    }
    else {
      if(length(main) != length(which))
        stop("number of plot titles is not correct")
      main2 <- character(6)
      main2[show] <- main
      main <- main2
    }
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        bvcpp(x, mar = 1, ci = ci, cilwd = cilwd, main = main[1],
              xlim = c(0,1), ylim = c(0,1), ...)
    }
    if (show[2]) {
        bvcpp(x, mar = 2, ci = ci, cilwd = cilwd, main = main[2],
              xlim = c(0,1), ylim = c(0,1), ...)
    }
    if (show[3]) {
        bvdens(x, grid = grid, legend = legend, main = main[3], ...)
    }
    if (show[4]) {
        bvdp(x, nplty = nplty, blty = blty, method = method,
             convex = convex, rev = rev, main = main[4], ...)
    }
    if (show[5]) {
        bvqc(x, p = p, mint = mint, legend = legend, main = main[5], ...)
    }
    if (show[6]) {
        bvh(x, half = half, main = main[6], ...)
    }
    invisible(x)
}

"bvcpp" <- function (x, ...) UseMethod("bvcpp")
"bvdens" <- function (x, ...) UseMethod("bvdens")
"bvdp" <- function (x, ...) UseMethod("bvdp")
"bvqc" <- function (x, ...) UseMethod("bvqc")
"bvh" <- function (x, ...) UseMethod("bvh")

"bvcpp.bvevd" <-  function(x, mar = 2, ci = TRUE, cilwd = 1, main = "Conditional Probability Plot", xlab = "Empirical", ylab = "Model", ...)
{ 
    data <- x$tdata
    mle.m1 <- x$param[c("loc1","scale1","shape1")]
    mle.m2 <- x$param[c("loc2","scale2","shape2")]
    data[,1:2] <- exp(-mtransform(data[,1:2], list(mle.m1, mle.m2)))
    narow <- is.na(data[,1]) | is.na(data[,2])
    data <- data[!narow,, drop=FALSE]
    n <- nrow(data)
    ppx <- ppoints(n)
    if(x$model %in% c("log","hr","neglog")) {
      probs <- ccbvevd(data, mar = mar, dep = x$param["dep"],
                    model = x$model)}
    if(x$model  %in% c("alog","aneglog"))
      probs <- ccbvevd(data, mar = mar, dep = x$param["dep"],
                    asy = x$param[c("asy1","asy2")], model = x$model)
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
      probs <- ccbvevd(data, mar = mar, alpha = x$param["alpha"],
                    beta = x$param["beta"], model = x$model)
    probs <- sort(probs)
    if(!ci) {
        plot(ppx, probs, main = main, xlab = xlab, ylab = ylab, ...)
        abline(0, 1)
    }
    else {
        samp <- runif(n*99)
        samp <- matrix(samp, n, 99)
        samp <- apply(samp, 2, sort)
        samp <- apply(samp, 1, sort)
        env <- t(samp[c(3,97),])
        matplot(ppx, cbind(probs, env), main = main, xlab = xlab,
                ylab = ylab, type = "pnn", pch = 4, ...)
        xyuser <- par("usr")
        smidge <- min(diff(c(xyuser[1], ppx, xyuser[2])))/2
        smidge <- max(smidge, (xyuser[2] - xyuser[1])/200)
        segments(ppx-smidge, env[,1], ppx+smidge, env[,1], lwd = cilwd)
        segments(ppx-smidge, env[,2], ppx+smidge, env[,2], lwd = cilwd)
        abline(0, 1)
    }
    invisible(list(x = ppx, y = probs))
}

"bvdens.bvevd" <-  function(x, grid = 50, legend = TRUE, pch = 1, main = "Density Plot", xlab = colnames(x$data)[1], ylab = colnames(x$data)[2], ...)
{
    xlimit <- range(x$tdata[,1], na.rm = TRUE)
    ylimit <- range(x$tdata[,2], na.rm = TRUE)
    xlimit[1] <- xlimit[1] - diff(xlimit) * 0.1
    xlimit[2] <- xlimit[2] + diff(xlimit) * 0.1
    ylimit[1] <- ylimit[1] - diff(ylimit) * 0.1
    ylimit[2] <- ylimit[2] + diff(ylimit) * 0.1
    xvec <- seq(xlimit[1], xlimit[2], length = grid)
    yvec <- seq(ylimit[1], ylimit[2], length = grid)
    xyvals <- expand.grid(xvec, yvec)
    mar1 <- x$param[c("loc1","scale1","shape1")]
    mar2 <- x$param[c("loc2","scale2","shape2")]
    if(x$model %in% c("log","hr","neglog"))
        dfunargs <- list(dep = x$param["dep"], mar1 = mar1, mar2 = mar2)
    if(x$model  %in% c("alog","aneglog"))
        dfunargs <- list(dep = x$param["dep"],
            asy = x$param[c("asy1","asy2")], mar1 = mar1, mar2 = mar2)
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
        dfunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"],
            mar1 = mar1, mar2 = mar2)
    dfunargs <- c(list(x = xyvals, model = x$model), dfunargs)
    dens <- do.call("dbvevd", dfunargs)
    dens <- matrix(dens, nrow = grid, ncol = grid)
    contour(xvec, yvec, dens, main = main, xlab = xlab, ylab = ylab, ...)
    data <- x$tdata
    if(ncol(data) == 2) points(data, pch = pch)
    if(ncol(data) == 3) {
      si <- data[,3] ; data <- data[,1:2]
      points(data[is.na(si),], pch = 4)
      points(data[si & !is.na(si),], pch = 16)
      points(data[!si & !is.na(si),], pch = 1)
      legwrd <- c("True","False","Unknown") ; legpch <- c(16,1,4)
      if(!any(is.na(si))) {legwrd <- legwrd[1:2] ; legpch <- legpch[1:2]}
      if(legend) legend(xlimit[1], ylimit[2], legwrd, pch = legpch)
    }
    invisible(list(x = xyvals, y = dens))
}

"bvdp.bvevd" <- function(x, method = "cfg", convex = FALSE, rev = FALSE, add = FALSE, lty = 1, nplty = 2, blty = 3, main = "Dependence Function", xlab = "t", ylab = "A(t)", ...)
{
    if(ncol(x$data) == 3) nplty <- 0
    abvnonpar(data = x$data[,1:2], nsloc1 = x$nsloc1, nsloc2 = x$nsloc2,
              epmar = FALSE, method = method, convex = convex, rev = rev,
              plot = TRUE, lty = nplty, blty = blty, main = main, xlab = xlab,
              ylab = ylab, add = add, ...)
    if(x$model %in% c("log","hr","neglog"))
        afunargs <- list(dep = x$param["dep"])
    if(x$model  %in% c("alog","aneglog"))
        afunargs <- list(dep = x$param["dep"], asy = x$param[c("asy1","asy2")])
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
        afunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"])
    afunargs <- c(list(rev = rev, add = TRUE, lty = lty, model = x$model),
                  afunargs)
    do.call("abvevd", afunargs)
    invisible(x)
}

"bvh.bvevd" <- function(x, half = FALSE, add = FALSE, lty = 1, main = "Spectral Density", xlab = "t", ylab = "h(t)", ...)
{
    if(x$model %in% c("log","hr","neglog"))
        afunargs <- list(dep = x$param["dep"])
    if(x$model  %in% c("alog","aneglog"))
        afunargs <- list(dep = x$param["dep"], asy = x$param[c("asy1","asy2")])
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
        afunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"])
    afunargs <- c(list(half = half, add = add, plot = TRUE, lty = lty, main = main,
       xlab = xlab, ylab = ylab, model = x$model), afunargs)
    do.call("hbvevd", afunargs)
    invisible(x)
}

"bvqc.bvevd"<- 
function(x, p = seq(0.75, 0.95, 0.05), mint = 1, add = FALSE,
   legend = TRUE, lty = 1, lwd = 1, col = 1, xlim = range(x$tdata[,1],
   na.rm = TRUE), ylim = range(x$tdata[,2], na.rm = TRUE), xlab =
   colnames(x$data)[1], ylab = colnames(x$data)[2], ...)
{
    if(mode(p) != "numeric" || p <= 0 || p >= 1)
      stop("`p' must be a vector of probabilities")
    nom <- 100
    om <- seq(0, 1, length = nom)
    # Calculate A(t)
    if(x$model %in% c("log","hr","neglog"))
        afunargs <- list(dep = x$param["dep"])
    if(x$model  %in% c("alog","aneglog"))
        afunargs <- list(dep = x$param["dep"], asy = x$param[c("asy1","asy2")])
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
        afunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"])
    afunargs <- c(list(x = om, plot = FALSE, model = x$model), afunargs)
    aom <- do.call("abvevd", afunargs)
    # End Calculate A(t)
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
      qct[[i]] <- -cbind(om/aom * log(p[i]), (1-om)/aom * log(p[i]))
      mar1 <- x$param[c("loc1","scale1","shape1")]
      mar2 <- x$param[c("loc2","scale2","shape2")]
      qct[[i]] <- mtransform(qct[[i]], list(mar1, mar2), inv = TRUE)    
      qct[[i]][1,1] <- 1.5 * xlim[2]
      qct[[i]][nom,2] <- 1.5 * ylim[2]
    }
    
    if(!add) {
      if(is.null(xlab)) xlab <- ""
      if(is.null(ylab)) ylab <- ""
      if(ncol(x$tdata) == 2) {
        plot(x$tdata[,1:2], xlab = xlab, ylab = ylab, xlim = xlim,
        ylim = ylim, ...)
      }
      if(ncol(x$tdata) == 3) {
        plot(x$tdata[,1:2], xlab = xlab, ylab = ylab, xlim = xlim,
        ylim = ylim, type = "n", ...)
        si <- x$tdata[,3] ; data <- x$tdata[,1:2]
        points(data[is.na(si),], pch = 4)
        points(data[si & !is.na(si),], pch = 16)
        points(data[!si & !is.na(si),], pch = 1)
        legwrd <- c("True","False","Unknown") ; legpch <- c(16,1,4)
        if(!any(is.na(si))) {legwrd <- legwrd[1:2] ; legpch <- legpch[1:2]}
        if(legend) legend(xlim[1], ylim[2], legwrd, pch = legpch)
      }
      for(i in 1:np) lines(qct[[i]], lty = lty, lwd = lwd, col = col)
    }
    else {
      for(i in 1:np) lines(qct[[i]], lty = lty, lwd = lwd, col = col)
    }
    return(invisible(qct))
}

### Bivariate POT Models ###

"plot.bvpot" <-  function(x, mar = 0, which = 1:4, main,
     ask = nb.fig < length(which) && dev.interactive(), grid = 50,
     above = FALSE, levels = NULL, tlty = 1, blty = 3, rev = FALSE,
     p = seq(0.75, 0.95, 0.05), half = FALSE, ...) 
{
    if (!inherits(x, "bvpot")) 
        stop("Use only with `bvpot' objects")
    nb.fig <- prod(par("mfcol"))
    if(mar == 1 || mar == 2) {
        indx <- paste(c("scale","shape"), as.character(mar), sep="")
        param <- x$param[indx]
        names(param) <- c("scale","shape")
        mdata <- x$data[, mar]
        mth <- x$threshold[mar]
        mexceed <- mdata[mdata > mth & !is.na(mdata)]       
        pot.mar <- structure(list(param = param, data = mdata,
           threshold = mth, exceedances = mexceed, nhigh = length(mexceed),
           npp = length(mdata), scale = param["scale"]),
           class = c("pot", "uvevd", "evd"))
        if(missing(which)) which <- 1:4
        plot(pot.mar, which = which, main = main, ask = ask, ...)
        return(invisible(x))
    }
    if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    if(missing(main)) {
      main <- c("Density Plot", "Dependence Function", "Quantile Curves",
                "Spectral Density")
    }
    else {
      if(length(main) != length(which))
        stop("number of plot titles is not correct")
      main2 <- character(4)
      main2[show] <- main
      main <- main2
    }
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        bvdens(x, grid = grid, above = above, levels = levels,
               tlty = tlty, main = main[1], ...)
    }
    if (show[2]) {
        bvdp(x, blty = blty, rev = rev, main = main[2], ...)
    }
    if (show[3]) {
        bvqc(x, p = p, above = above, tlty = tlty, main = main[3], ...)
    }
    if (show[4]) {
        bvh(x, half = half, main = main[4], ...)
    }
    invisible(x)
}

"bvdens.bvpot" <-  function(x, grid = 50, above = FALSE, tlty = 1, levels = NULL, main = "Density Plot", pch = 1, xlab = colnames(x$data)[1], ylab = colnames(x$data)[2], xlim, ylim, ...)
{
    xlimit <- range(x$data[,1], na.rm = TRUE)
    ylimit <- range(x$data[,2], na.rm = TRUE)
    xlimit[1] <- xlimit[1] - diff(xlimit) * 0.1
    xlimit[2] <- xlimit[2] + diff(xlimit) * 0.1
    ylimit[1] <- ylimit[1] - diff(ylimit) * 0.1
    ylimit[2] <- ylimit[2] + diff(ylimit) * 0.1
    if(missing(xlim)) xlim <- xlimit
    if(missing(ylim)) ylim <- ylimit
    u1 <- x$threshold[1]
    u2 <- x$threshold[2]
    if((xlimit[2] <= u1) || (ylimit[2] <= u2))
      stop("x and y limits must contain thresholds")
    xvec <- seq(u1, xlimit[2], length = grid)
    yvec <- seq(u2, ylimit[2], length = grid)
    xyvals <- txyvals <- fxyvals <- expand.grid(xvec, yvec)
    mar1 <- x$param[c("scale1","shape1")]
    mar2 <- x$param[c("scale2","shape2")]
    # Transform exceedance grid to frechet margins
    txyvals[,1] <- mtransform(xyvals[,1], c(u1, mar1))
    txyvals[,2] <- mtransform(xyvals[,2], c(u2, mar2))
    lambda <- x$nat[1:2]/(nrow(x$data) + 1)
    fxyvals[,1] <- -1/log(1 - lambda[1] * txyvals[,1])
    fxyvals[,2] <- -1/log(1 - lambda[2] * txyvals[,2])
    # End transform
    if(x$model %in% c("log","hr","neglog"))
        dfunargs <- list(dep = x$param["dep"], mar1 = c(1,1,1), mar2 = c(1,1,1))
    if(x$model  %in% c("alog","aneglog"))
        dfunargs <- list(dep = x$param["dep"],
            asy = x$param[c("asy1","asy2")], mar1 = c(1,1,1), mar2 = c(1,1,1))
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
        dfunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"],
            mar1 = c(1,1,1), mar2 = c(1,1,1))
    dfunargs <- c(list(x = fxyvals, model = x$model), dfunargs)
    # Jacobian terms
    txyvals[,1] <- fxyvals[,1]^2 * txyvals[,1]^(1 + mar1[2]) /
      (1 - lambda[1] * txyvals[,1])
    txyvals[,1] <- lambda[1] * txyvals[,1] / mar1[1]
    txyvals[,2] <- fxyvals[,2]^2 * txyvals[,2]^(1 + mar2[2]) /
      (1 - lambda[2] * txyvals[,2])
    txyvals[,2] <- lambda[2] * txyvals[,2] / mar2[1]
    # End jacobian terms 
    dens <- do.call("dbvevd", dfunargs)
    dens <- dens * txyvals[,1] * txyvals[,2]
    dens <- matrix(dens, nrow = grid, ncol = grid)
    if(is.null(levels)) {
      levels <- seq(10, 40, length = 4)
      levels <- dens[cbind(levels, levels)]
      levels <- signif(levels, 1)
    }
    contour(xvec, yvec, dens, main = main, xlab = xlab, ylab = ylab,
            xlim = xlim, ylim = ylim, levels = levels, ...)
    abline(v = u1, lty = tlty); abline(h = u2, lty = tlty)
    data <- x$data
    if(above) {
      above <- (data[,1] > u1) & (data[,2] > u2)
      data <- data[above,]
    }
    points(data, pch = pch)
    invisible(list(x = xyvals, y = dens))
}

"bvdp.bvpot" <- function(x, rev = FALSE, add = FALSE, lty = 1, blty = 3, main = "Dependence Function", xlab = "t", ylab = "A(t)", ...)
{
    if(x$model %in% c("log","hr","neglog"))
        afunargs <- list(dep = x$param["dep"])
    if(x$model  %in% c("alog","aneglog"))
        afunargs <- list(dep = x$param["dep"], asy = x$param[c("asy1","asy2")])
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
        afunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"])
    afunargs <- c(list(rev = rev, add = add, plot = TRUE, lty = lty, blty = blty,
      main = main, xlab = xlab, ylab = ylab, model = x$model), afunargs)
    do.call("abvevd", afunargs)
    invisible(x)
}

"bvqc.bvpot"<- 
function(x, p = seq(0.75, 0.95, 0.05), above = FALSE, tlty = 1,
   add = FALSE, lty = 1, lwd = 1, col = 1, xlim =
   range(x$data[,1], na.rm = TRUE), ylim =
   range(x$data[,2], na.rm = TRUE), xlab = colnames(x$data)[1],
   ylab = colnames(x$data)[2], ...)
{
    if(mode(p) != "numeric" || p <= 0 || p >= 1)
      stop("`p' must be a vector of probabilities")
    nom <- 100
    om <- seq(0, 1, length = nom)
    # Calculate A(t)
    if(x$model %in% c("log","hr","neglog"))
        afunargs <- list(dep = x$param["dep"])
    if(x$model  %in% c("alog","aneglog"))
        afunargs <- list(dep = x$param["dep"], asy = x$param[c("asy1","asy2")])
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
        afunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"])
    afunargs <- c(list(x = om, plot = FALSE, model = x$model), afunargs)
    aom <- do.call("abvevd", afunargs)
    # End Calculate A(t)
    np <- length(p)
    qct <- list()
    if(add) {
      xlim <- par("usr")[1:2]
      ylim <- par("usr")[1:2]
      if(par("xlog")) xlim <- 10^xlim
      if(par("ylog")) ylim <- 10^ylim
    }
    u1 <- x$threshold[1]
    u2 <- x$threshold[2]
    lambda <- x$nat[1:2]/(nrow(x$data) + 1)
    for(i in 1:np) {
      qct[[i]] <- cbind((1 - p[i]^(om/aom))/lambda[1],
         (1 - p[i]^((1-om)/aom))/lambda[2])
      mar1 <- c(u1, x$param[c("scale1","shape1")])
      mar2 <- c(u2, x$param[c("scale2","shape2")])
      qct[[i]] <- mtransform(qct[[i]], list(mar1, mar2), inv = TRUE)    
      qct[[i]][1,1] <- 1.5 * xlim[2]
      qct[[i]][nom,2] <- 1.5 * ylim[2]
    }
    data <- x$data
    if(above) {
      above <- (data[,1] > u1) & (data[,2] > u2)
      data <- data[above,]
    }
    if(!add) {
      if(is.null(xlab)) xlab <- ""
      if(is.null(ylab)) ylab <- ""
      plot(data, xlab = xlab, ylab = ylab, xlim = xlim,
           ylim = ylim, ...)
      abline(v = u1, lty = tlty); abline(h = u2, lty = tlty)
      for(i in 1:np) lines(qct[[i]], lty = lty, lwd = lwd, col = col)
    }
    else {
      for(i in 1:np) lines(qct[[i]], lty = lty, lwd = lwd, col = col)
    }
    return(invisible(qct))
}

"bvh.bvpot" <- function(x, half = FALSE, add = FALSE, lty = 1, main = "Spectral Density", xlab = "t", ylab = "h(t)", ...)
{
    if(x$model %in% c("log","hr","neglog"))
        afunargs <- list(dep = x$param["dep"])
    if(x$model  %in% c("alog","aneglog"))
        afunargs <- list(dep = x$param["dep"], asy = x$param[c("asy1","asy2")])
    if(x$model  %in% c("bilog","negbilog","ct","amix"))
        afunargs <- list(alpha = x$param["alpha"], beta = x$param["beta"])
    afunargs <- c(list(half = half, add = add, plot = TRUE, lty = lty,
       main = main, xlab = xlab, ylab = ylab, model = x$model), afunargs)
    do.call("hbvevd", afunargs)
    invisible(x)
}

### Documented Ancillary Functions ###

"mtransform"<-
function(x, p, inv = FALSE, drp = FALSE)
{
    if(is.list(p)) {
      if(is.null(dim(x)) && length(x) != length(p))
        stop(paste("`p' must have", length(x), "elements"))
      if(!is.null(dim(x)) && ncol(x) != length(p))
        stop(paste("`p' must have", ncol(x), "elements"))
      if(is.null(dim(x))) dim(x) <- c(1, length(p))
      for(i in 1:length(p))
        x[,i] <- Recall(x[,i], p[[i]], inv = inv)
      if(ncol(x) == 1 || (nrow(x) == 1 && drp)) x <- drop(x)
      return(x)
    }
    if(is.null(dim(x))) dim(x) <- c(length(x), 1)
    p <- matrix(t(p), nrow = nrow(x), ncol = 3, byrow = TRUE)
    if(min(p[,2]) <= 0) stop("invalid marginal scale")
    gumind <- (p[,3] == 0)
    nzshapes <- p[!gumind,3]
    if(!inv) {
        x <- (x - p[,1])/p[,2]
        x[gumind, ] <- exp(-x[gumind, ])
        if(any(!gumind))
        x[!gumind, ] <- pmax(1 + nzshapes*x[!gumind, ], 0)^(-1/nzshapes)
    }
    else {
        x[gumind, ] <- p[gumind,1] - p[gumind,2] * log(x[gumind, ])
        x[!gumind, ] <- p[!gumind,1] + p[!gumind,2] *
          (x[!gumind, ]^(-nzshapes) - 1)/nzshapes
    }
    if(ncol(x) == 1 || (nrow(x) == 1 && drp)) x <- drop(x)
    x
}

"ccbvevd" <- function(x, mar = 2, dep, asy = c(1, 1), alpha, beta, model =
  c("log", "alog", "hr", "neglog", "aneglog", "bilog", "negbilog", "ct",
  "amix"), lower.tail = TRUE)
{
  if(min(x[,1:2]) <= 0 || max(x[,1:2]) >= 1)
    stop("x must contain values in (0,1)")
  model <- match.arg(model)
  m1 <- c("bilog", "negbilog", "ct", "amix")
  m2 <- c(m1, "log", "hr", "neglog")
  m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  if((model %in% m1) && !missing(dep))
    warning("ignoring `dep' argument")
  if((model %in% m2) && !missing(asy)) 
    warning("ignoring `asy' argument")
  if((model %in% m3) && !missing(alpha))
    warning("ignoring `alpha' argument")
  if((model %in% m3) && !missing(beta)) 
    warning("ignoring `beta' argument")
  if(model %in% m1) dep <- 1
  if(model %in% m3) alpha <- beta <- 1
  imodel <- match(model, c("log","alog","hr","neglog","aneglog",
    "bilog","negbilog","ct","amix"))
  n <- nrow(x)
  
  if(ncol(x) == 2) {
    ccop <- .C("ccop", as.double(x[,1]), as.double(x[,2]), as.integer(mar),
      as.double(dep), as.double(asy[1]), as.double(asy[2]), as.double(alpha),
      as.double(beta), as.integer(n), as.integer(imodel), ccop = double(n),
      PACKAGE = "evd")$ccop
  }
  
  if(ncol(x) == 3) {

    "dbvevd.case" <- function(x1, x2, case, mar, dep, asy, alpha, beta)
    {
      nlbvfn <- paste("nlbv", model, sep = "")
      n <- max(length(x1), length(x2))
      x1 <- rep(-1/log(x1), length = n)
      x2 <- rep(-1/log(x2), length = n)
      case <- rep(case, length = n)
      mpar <- as.double(1)
      split <- as.integer(1)
      if(mar == 1) {
        tmp <- x1; x1 <- x2; x2 <- tmp
        if(model %in% c("alog","aneglog")) asy <- rev(asy)
        if(model %in% c("bilog","negbilog","ct"))
          { tmp <- alpha; alpha <- beta; beta <- tmp }
        if(model == "amix")
          { alpha <- alpha + 3*beta; beta <- -beta }
      }
      if(model %in% c("bilog","negbilog","ct","amix"))
        nl <- .C(nlbvfn,
          as.double(x1), as.double(x2), n, case,
          as.double(alpha), as.double(beta), rep(mpar,n), mpar, mpar,
          rep(mpar,n), mpar, mpar, split, dns = double(n),
          PACKAGE = "evd")$dns
      if(model %in% c("log","hr","neglog"))
        nl <- .C(nlbvfn,
          as.double(x1), as.double(x2), n, case,
          as.double(dep), rep(mpar,n), mpar, mpar, rep(mpar,n), mpar,
          mpar, split, dns = double(n), PACKAGE = "evd")$dns
      if(model %in% c("alog","aneglog"))
        nl <- .C(nlbvfn,
          as.double(x1), as.double(x2), n, case,
          as.double(dep), as.double(asy[1]), as.double(asy[2]), rep(mpar,n),
          mpar, mpar, rep(mpar,n), mpar, mpar, split, dns = double(n),
          PACKAGE = "evd")$dns
      jac.alt <- 1/x1 + 1/x2 + 2*log(x1 * x2)
      exp(jac.alt - nl)
    }

    ccop <- numeric(n)
    case <- as.integer(x[,3])
    eps <- .Machine$double.eps^0.5
    if(mar == 2) { fm <- x[,1] ; cm <- x[,2] }
    if(mar == 1) { fm <- x[,2] ; cm <- x[,1] }
    for(i in 1:n) {
      if(is.na(case[i])) {
        ccop[i] <- .C("ccop", as.double(x[i,1]), as.double(x[i,2]),
          as.integer(mar), as.double(dep), as.double(asy[1]),
          as.double(asy[2]), as.double(alpha), as.double(beta),
          as.integer(1), as.integer(imodel), ccop = double(1),
          PACKAGE = "evd")$ccop
      }
      else {
        den <- integrate("dbvevd.case", eps, 1-eps, x2 = cm[i], case =
          case[i], mar=mar, dep=dep, asy=asy, alpha=alpha, beta=beta)$value
        num <- integrate("dbvevd.case", eps, fm[i], x2 = cm[i], case =
          case[i], mar=mar, dep=dep, asy=asy, alpha=alpha, beta=beta)$value
        ccop[i] <- num/den
      }
    }
  }
  if(!lower.tail) ccop <- 1 - ccop
  ccop
}





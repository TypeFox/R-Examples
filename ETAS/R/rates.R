
rates <- function(fit, dimyx=NULL, method="zhuang", plot.it=TRUE)
{
  spatstat::verifyclass(fit, "etas")
  rates.inter(fit$param, fit$object, fit$bwd, dimyx=dimyx,
              method=method, plot.it=plot.it)
}

rates.inter <- function(theta, object, bwd, dimyx=NULL,
                  method="zhuang", plot.it=TRUE)
{
  xx <- object$longlat$long
  yy <- object$longlat$lat
  win <- spatstat::owin(poly=list(x=object$region.poly$long,
                                  y=object$region.poly$lat))
  tperiod <- object$rtperiod
  revents <- object$revents
  bk <- revents[, 6]
  pb <- revents[, 7]
  lam <- revents[, 8]

  switch(method, spatstat={
    Xg <- spatstat::ppp(xx, yy, window=spatstat::owin(range(xx), range(yy)))
    total.im <- density(Xg, weights=rep(1/diff(tperiod), Xg$n), dimyx)[win]
    spatstat::marks(Xg) <- bk
    bkgd.im <-  spatstat::smooth.ppp(Xg, weights=bk/diff(tperiod))[win]
    spatstat::marks(Xg) <- lam - bk
    clust.im <- spatstat::smooth.ppp(Xg,
                                     weights=rep(1/diff(tperiod), Xg$n))[win]
  }, zhuang={
    if (is.null(dimyx))
      dimyx <- c(128, 128)
    if (!is.numeric(dimyx) || length(dimyx) != 2)
      stop(paste(sQuote(dimyx), "must be a numeric vector of length 2."))
    gr <- spatstat::gridcenters(spatstat::boundingbox(win), dimyx[2], dimyx[1])
    gx <- gr$x
    gy <- gr$y
    out1 <- out2 <- out3 <- numeric(length(gx))
    for (i in 1:length(gx))
    {
      r2 <- (xx - gx[i])^2 + (yy - gy[i])^2
      s1 <- exp(-r2/(2 * bwd^2)) / (2 * pi * bwd^2)
      s2 <- pb *  s1
      s1 <- sum(s1)/diff(tperiod)
      s2 <- sum(s2)/diff(tperiod)
      out1[i] <- s1
      out2[i] <- s2
      out3[i] <- 1 - s2 / s1
    }
    total.im <- spatstat::as.im(list(x=unique(gx), y=unique(gy),
                                     z=matrix(out1, nrow=dimyx[2], ncol=dimyx[1])))
    bkgd.im <-  spatstat::as.im(list(x=unique(gx), y=unique(gy),
                                     z=matrix(out2, nrow=dimyx[2], ncol=dimyx[1])))
    clust.im <- spatstat::as.im(list(x=unique(gx), y=unique(gy),
                                     z=matrix(out3, nrow=dimyx[2], ncol=dimyx[1])))
  })


  out <- list(total=total.im, background=bkgd.im, clustring=clust.im)
  if (plot.it)
  {
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(1, 2), mar=c(1.5, 1.5, 2, 2))
    rp <- list(x=object$region.poly$long, y=object$region.poly$lat)
    plot(bkgd.im[win, drop=FALSE], main="background seismicity rate", axes=TRUE, ribsep=0.1)
    map('world', add=TRUE, col="grey50")
    plot(clust.im[win, drop=FALSE], main="clustering coefficient", axes=TRUE, ribsep=0.1)
    map('world', add=TRUE, col="grey50")
    par(oldpar)
    invisible(out)
  }
  else
    return(out)
}

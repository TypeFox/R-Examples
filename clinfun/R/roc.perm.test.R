roc.perm.test <- function(marker, status, marker2=NULL, group=NULL, nperm=2500, mp=NULL) {
  if (any(!is.finite(marker))) stop("Marker values should be finite")
  if (any(!is.finite(status))) stop("All status should be finite")
  if (sum(sapply(list(marker2, group), is.null)) != 1) stop("exactly one of marker2 or group must be specified")
  if (!is.null(marker2)) {
    if (any(!is.finite(marker2))) stop("Marker2 values should be finite")
    out <- compareROC.paired(marker,marker2,status,nperm)
  } else {
    ug <- unique(group)
    if (length(ug) !=2) stop("group should have exactly 2 categories")
    x <- marker[group == ug[1]]
    dx <- status[group == ug[1]]
    y <- marker[group == ug[2]]
    dy <- status[group == ug[2]]
    out <- compareROC.unpaired(x, dx, y, dy, nperm, mp)
  }
  out
}

print.roc.perm.test <- function(x, ...) {
  if (!inherits(x, 'roc.perm.test')) stop("Object not of class roc.perm.test")
  cat("Statistic: ", x$ostat, "; p-value =", x$p.value, "; based on", length(x$pstat), "permutations \n")
}

plot.roc.perm.test <- function(x, ...) {
  if (!inherits(x, 'roc.perm.test')) stop("Object not of class roc.perm.test")
  oo <- density(c(x$ostat,x$pstat))
  plot(oo, xlab="Test statistic", main="", ...)
  oy <- max(oo$y)
  arrows(x$ostat, 0.15*oy, x$ostat, 0, length=0.12, lwd=1.5)
}

compareROC.paired <- function(x,y,d,nperm=2500){
  n <- length(x)
  rx <- rank(x)
  ry <- rank(y)
  dx <- d[order(rx)]
  dy <- d[order(ry)]
  prx <- pry <- pdx <- pdy <- uu <- rep(0,n)
  ostat <-  2*sum(abs(cumsum(dx) - cumsum(dy)))/n^2
  pstat <- sapply(1:nperm,function(i,n,rx,ry,d,prx,pry,pdx,pdy,uu) {
    prx <- rx + runif(n) - 0.5
    pry <- ry + runif(n) - 0.5
    uu <- 1*(runif(n) < 0.5)
    pdx <- d[order(prx*uu + pry*(1-uu))]
    pdy <- d[order(prx*(1-uu) + pry*uu)]
    2*sum(abs(cumsum(pdx) - cumsum(pdy)))/n^2
  },n,rx,ry,d,prx,pry,pdx,pdy,uu)
  out <- list()
  out$ostat <- ostat
  out$pstat <- pstat
  out$p.value <- sum(pstat>=ostat)/nperm
  class(out) <- "roc.perm.test"
  out
}

compareROC.unpaired <- function(x, dx, y, dy, nperm=2500, mp=NULL){
  nx <- length(x)
  nx1 <- sum(dx)
  nx0 <- nx - nx1
  ny <- length(y)
  ny1 <- sum(dy)
  ny0 <- ny - ny1
  if (is.null(mp) | missing(mp)) mp <- (nx1 + ny1)/(nx + ny)
  one <- rep(1, nx+ny)
  dx <- dxp <- dx[order(x)]
  dy <- dyp <- dy[order(y)]
  xx <- cumsum((1-mp)*(1-dx)/nx0 + mp*dx/nx1)
  yy <- cumsum((1-mp)*(1-dy)/ny0 + mp*dy/ny1)
  xy <- c(xx,yy)
  gg <- rep(0:1, c(nx, ny))
  oxy <- order(c(xx,yy))
  xy <-  pxy <- xy[oxy]
  gxy <- gxyp <- gg[oxy]
  dxy <- c(dx,dy)[oxy]
  ex <- pex <- (1-mp)*(one-cumsum((1-gxy)*(1-dxy)/nx0)) + mp*cumsum((1-gxy)*dxy/nx1)
  ey <- pey <- (1-mp)*(one-cumsum(gxy*(1-dxy)/ny0)) + mp*cumsum(gxy*dxy/ny1)
  ostat <- sum(diff(c(0,xy))*abs(ex-ey))
  pstat <- sapply(1:nperm, function(i, pg, pd, one, dxp, dyp, nx0, nx1, ny0, ny1, mp, pxy, pex, pey) {
    pg[pd==0] <- sample(pg[pd==0])
    pg[pd==1] <- sample(pg[pd==1])
    dxp <- pd[pg==0]
    dyp <- pd[pg==1]
    pxy <- sort(c(cumsum((1-mp)*(1-dxp)/nx0 + mp*dxp/nx1), cumsum((1-mp)*(1-dyp)/ny0 + mp*dyp/ny1)))
    pex <- (1-mp)*(one-cumsum((1-pg)*(1-pd)/nx0)) + mp*cumsum((1-pg)*pd/nx1)
    pey <- (1-mp)*(one-cumsum(pg*(1-pd)/ny0)) + mp*cumsum(pg*pd/ny1)
    sum(diff(c(0,pxy))*abs(pex-pey))
  }, gxy, dxy, one, dxp, dyp, nx0, nx1, ny0, ny1, mp, pxy, pex, pey)
  out <- list()
  out$ostat <- ostat
  out$pstat <- pstat
  out$p.value <- sum(pstat>=ostat)/nperm
  class(out) <- "roc.perm.test"
  out
}

alternativesenv <-
function (dat, g1, group1 = "Species", ev, vars, world, xname = "", 
          yname = "", rst, locality = "", ext = c(-180, 180, -60, 90)) 
{
  fg <- dev.list()
  fig1 <- fg[1]
  fig2 <- fg[2]
  dev.off(which = fg[3])
  length(dat)
  length(world)
  length(rst)
  fieldsmissing(dat, fields = c("ID", "x", "y", "x_original", 
                                "y_original", "Correction", "Exclude","Reason"))
  cn <- names(dat)
  f3 <- match(group1, cn)
  ff <- which(dat[, f3] == g1 & dat$Exclude == 0)
  x1 <- coord2numeric(dat$x[ff])
  y1 <- coord2numeric(dat$y[ff])
  nas <- cbind(is.na(x1) * 1, is.na(y1) * 1)
  f <- which(nas >= 1)
  ff2 <- ff[f]
  dat$Exclude[ff2] <- 1
  ff <- which(dat[, f3] == g1 & dat$Exclude == 0)
  x <- coord2numeric(dat$x[ff])
  y <- coord2numeric(dat$y[ff])
  datid <- dat$ID[ff]
  xy1 <- data.frame(x = x, y = y)
  xy <- SpatialPoints(xy1)
  ex <- getextent(x, y, ext)
  xlm <- ex$xlm
  ylm <- ex$ylm
  #bringToTop(dev.set(fig1), stay = T)
  dev.set(fig1)
  plot(world, border = "gray", xlim = xlm, ylim = ylm, axes=T)
  dat$Correction <- as.character(dat$Correction)
  vals <- extract(rst, xy)
  f2 <- which(!is.na(vals))
  f3 <- which(is.na(vals))
  points(x[f3], y[f3], pch = 20, cex = 1, col = "red")
  points(x[f2], y[f2], pch = 20, cex = 1, col = "blue")
  x1 <- x
  y1 <- y
  fs <- identify(x, y, n = 1, plot = FALSE)
  if (length(fs) == 0) {
    graphics.off()
    stop("no match, function stopped")
  }
  idf <- datid[fs]
  ff <- which(dat$ID == idf)
  if (locality == "") {
  } else {
    f1 <- match(locality, cn)
    mtext(dat[ff, f1], side = 3, line = 2)
  }
  x <- x[fs]
  y <- y[fs]
  mtext(paste("x = ", x, "; y = ", y, sep = ""), side = 3, 
        line = 0.5, cex = 0.9)
  text(x, y, labels = idf, pos = 4, adj = 1, cex = 0.8)
  points(y, x, pch = 20, col = "purple")
  segments(x, y, y, x, lty = 2)
  points(x * -1, y, pch = 20, col = "purple")
  segments(x, y, x * -1, y, lty = 2)
  points(x, y * -1, pch = 20, col = "purple")
  segments(x, y, x, y * -1, lty = 2)
  points(x * -1, y * -1, pch = 20, col = "purple")
  segments(x, y, x * -1, y * -1, lty = 2)
  xs <- substddmm(x)
  ys <- substddmm(y)
  points(xs, ys, pch = 20, col = "purple")
  segments(x, y, xs, ys, lty = 2)
  signx <- ifelse(x < 0, -1, 1)
  signy <- ifelse(y < 0, -1, 1)
  xs2 <- abs(y) * signx
  ys2 <- abs(x) * signy
  points(xs2, ys2, pch = 20, col = "purple")
  segments(x, y, xs2, ys2, lty = 2)
  xn <- c(y, x * -1, x, x * -1, xs, xs2)
  yn <- c(x, y, y * -1, y * -1, ys, ys2)
  xnew <- c(xn, x, x1)
  ynew <- c(yn, y, y1)
  xn1 <- min(xn)
  xn1 <- ifelse(xn1<xlm[1], xn1, xlm[1])
  xn2 <- max(xn)
  xn2 <- ifelse(xn2>xlm[2], xn2, xlm[2])
  xlmnew <- c(xn1, xn2)
  yn1 <- min(yn)
  yn1 <- ifelse(yn1<ylm[1], yn1, ylm[1])
  yn2 <- max(yn)
  yn2 <- ifelse(yn2>ylm[2], yn2, ylm[2])
  ylmnew <- c(yn1, yn2)
  plot(world, border = "gray", xlim = xlmnew, ylim = ylmnew, axes=T)
  points(x[f3], y[f3], pch = 20, cex = 1, col = "red")
  points(x[f2], y[f2], pch = 20, cex = 1, col = "blue")
  if (locality == "") {
  } else {
    f1 <- match(locality, cn)
    mtext(dat[ff, f1], side = 3, line = 2)
  }
  mtext(paste("x = ", x, "; y = ", y, sep = ""), side = 3, 
        line = 0.5, cex = 0.9)
  text(x, y, labels = idf, pos = 4, adj = 1, cex = 0.8)
  points(y, x, pch = 20, col = "purple")
  segments(x, y, y, x, lty = 2)
  points(x * -1, y, pch = 20, col = "purple")
  segments(x, y, x * -1, y, lty = 2)
  points(x, y * -1, pch = 20, col = "purple")
  segments(x, y, x, y * -1, lty = 2)
  points(x * -1, y * -1, pch = 20, col = "purple")
  segments(x, y, x * -1, y * -1, lty = 2)
  points(xs, ys, pch = 20, col = "purple")
  segments(x, y, xs, ys, lty = 2)
  points(xs2, ys2, pch = 20, col = "purple")
  segments(x, y, xs2, ys2, lty = 2)
  text(xn, yn, labels = 1:6, pos = 4, adj = 0, cex = 0.7)
  xynew <- data.frame(x = xnew, y = ynew)
  ed <- extract(ev, xynew)
  ned <- names(as.data.frame(ed))
  nv <- length(vars)
  vrs <- rep(0, nv)
  for (i in 1:nv) {
    vrs[i] <- match(vars[i], ned)
  }
  env <- ed[, vrs]
  fsea <- which(is.na(env[1:6, 1]))
  points(xn[fsea], yn[fsea], pch = 20, col = "red")
  nw <- 1:6
  #bringToTop(dev.set(fig2), stay = T)
  dev.set(fig2)
  plot(env[, 1], env[, 2], col = "black", xlab = xname, ylab = yname)
  points(env[nw, 1], env[nw, 2], col = "blue", pch = 18)
  points(env[7, 1], env[7, 2], col = "red", pch = 18)
  text(env[nw, 1], env[nw, 2], labels = 1:6, adj = 1, pos = 4, 
       cex = 0.7)
  text(env[7, 1], env[7, 2], labels = idf, pos = 4, adj = 1, 
       cex = 0.8)
  #bringToTop(dev.set(fig1), stay = T)
  dev.set(fig1)
  f1 <- identify(xnew, ynew, n = 1, plot = F)
  if (length(f1) == 0) {
    graphics.off()
    stop("No match found, function Stopped")
  }
  x1 <- xnew[f1]
  y1 <- ynew[f1]
  points(x1, y1, col = "black")
  if (f1 == 7) {
    dat$Exclude[ff] <- 1
  }
  else {
    dc <- dat$Correction[ff]
    dc <- str_replace(dc, "[.]", as.character(f1))
    dat$Correction[ff] <- dc
    dat$x[ff] <- x1
    dat$y[ff] <- y1
    dat$x_original[ff] <- x
    dat$y_original[ff] <- y
    xyp <- paste("new x = ", dat$x[ff], "; y = ", dat$y[ff], 
                 sep = "")
    mtext(xyp, side = 1, line = 0, cex = 0.8)
  }
  return(dat)
}

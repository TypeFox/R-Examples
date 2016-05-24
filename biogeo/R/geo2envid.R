geo2envid <-
function (edat, g1, group1 = "Species", group2 = "", world, xc = "AP", 
          yc = "AMT", xname = "Annual Precipitation (mm)", yname = "Mean Annual Temperature", 
          showrecord = "", ext = c(-180, 180, -60, 90)) 
{
  fg <- dev.list()
  fig1 <- fg[1]
  fig2 <- fg[2]
  fig3 <- fg[3]
  cn <- names(edat)
  f1 <- match(xc, cn)
  f2 <- match(yc, cn)
  f3 <- match(group1, cn)
  ff <- which(edat[, f3] == g1 & edat$Exclude == 0)
  xc <- edat[ff, f1]
  yc <- edat[ff, f2]
  x <- coord2numeric(edat$x[ff])
  y <- coord2numeric(edat$y[ff])
  nas <- cbind(is.na(xc) * 1, is.na(yc) * 1, is.na(x) * 1, 
               is.na(y) * 1)
  f <- which(nas >= 1)
  ff2 <- ff[f]
  edat$Exclude[ff2] <- 2
  id <- edat$ID[ff]
  ex <- getextent(x, y, ext)
  fids <- {}
  xlm <- ex$xlm
  ylm <- ex$ylm
  #bringToTop(dev.set(fig1), stay = T)
  dev.set(fig1)
  par(mai = c(1.36, 0.2, 1.093333, 0.2))
  plot(world, border = "gray", xlim = xlm, ylim = ylm, axes=T)
  box(which = "plot", lty = "solid", col = "brown")
  points(x, y, xlab = "x-coordinate", ylab = "y-coordinate")
  if (nchar(showrecord) > 0) {
    fid <- which(edat$ID == showrecord)
    if (length(fid) > 0) {
      text(edat$x[fid], edat$y[fid], labels = showrecord, 
           adj = 1, pos = 4, cex = 0.7)
      points(edat$x[fid], edat$y[fid], col = "green")
    }
  }
  if (nchar(group2) > 0) {
    f4 <- match(group2, cn)
    if (is.na(f4)) 
      stop("field for group2 not found in dat")
    g2 <- unique(edat[, f4])
    ng2 <- length(g2)
    if (ng2 > 0) {
      cols <- rep(1, ng2)
      for (j in 1:ng2) {
        fg2 <- which(edat[ff, f4] == g2[j])
        points(x[fg2], y[fg2], col = j + 2)
        cols[j] <- j + 2
      }
      legend("bottomleft", legend = as.character(g2), col = cols, 
             pch = 18)
    }
  }
  mtext(g1, side = 3, line = 3, adj = 0.5)
  #bringToTop(dev.set(fig2), stay = T)
  dev.set(fig2)
  par(mai = c(1.36, 1.093333, 1.093333, 0.5))
  plot(xc, yc, xlab = xname, ylab = yname)
  box(which = "plot", lty = "solid", col = "blue")
  if (nchar(group2) > 0) {
    if (ng2 > 0) {
      for (j in 1:ng2) {
        fg2 <- which(edat[ff, f4] == g2[j])
        points(xc[fg2], yc[fg2], col = j + 2)
      }
    }
  }
  f <- 0
  while (f != 4) {
    #bringToTop(dev.set(fig3), stay = T)
    dev.set(fig3)
    par(mai = c(0, 0, 0, 0))
    plot(1, 1, type = "n", xlab = "", ylab = "", xaxt = "n", 
         yaxt = "n", xlim = c(0, 5), ylim = c(0, 4))
    x1 <- rep(0.5, 4)
    y1 <- (4:1) - 0.5
    points(x1, y1, pch = 18, col = c("brown", "blue", "blue", 
                                     "black"))
    text(x1, y1, labels = c("Query Geographical", "Query Environment", 
                            "Exclude Environment", "Exit"), adj = 0, pos = 4, 
         cex = 0.7)
    f <- identify(x1, y1, n = 1, plot = FALSE)
    if (f == 1) {
      #bringToTop(dev.set(fig1))
      dev.set(fig1)
      fi <- identify(x, y, n = 1, plot = FALSE)
      points(x[fi], y[fi], pch = 18, col = "red")
      fdat <- ff[fi]
      fids[fdat] <- 1
      text(x[fi], y[fi], labels = edat$ID[fdat], adj = 1, pos = 4, 
           cex = 0.7)
      #bringToTop(dev.set(fig2))
      dev.set(fig2)
      points(xc[fi], yc[fi], pch = 18, col = "red")
      text(xc[fi], yc[fi], labels = edat$ID[fdat], adj = 1, pos = 4, 
           cex = 0.7)
    }
    if (f == 2) {
      #bringToTop(dev.set(fig2))
      dev.set(fig2)
      fi <- identify(xc, yc, n = 1, plot = FALSE)
      points(xc[fi], yc[fi], pch = 18, col = "red")
      fdat <- ff[fi]
      fids[fdat] <- 1
      text(xc[fi], yc[fi], labels = edat$ID[fdat], adj = 1, pos = 4, 
           cex = 0.7)
      #bringToTop(dev.set(fig1))
      dev.set(fig1)
      points(x[fi], y[fi], pch = 18, col = "red")
      text(x[fi], y[fi], labels = edat$ID[fdat], adj = 1, pos = 4, 
           cex = 0.7)
    }
    if (f == 3) {
      #bringToTop(dev.set(fig2))
      dev.set(fig2)
      fi <- identify(xc, yc, n = 1, plot = FALSE)
      points(xc[fi], yc[fi], pch = 18, col = "black")
      fdat <- ff[fi]
      fids[fdat] <- 1
      text(xc[fi], yc[fi], labels = edat$ID[fdat], adj = 1, pos = 4, 
           cex = 0.7)
      edat$Exclude[fdat] <- 5
    }
    if (f == 4) {
      graphics.off()
    }
  }
  if (any(edat$Exclude == 2)) {
    f2 <- which(edat$Exclude == 2)
    edat$Exclude[f2] <- 0
  }
  if(length(fids[!is.na(fids)])>0){
    print(paste(edat$ID[!is.na(fids)],collapse='; '))
  }
  return(edat)
}

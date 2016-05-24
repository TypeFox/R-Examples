alternatives2 <-
function (dat, g1, group1 = "Species", group2 = "", world, rst, 
          locality = "", pos = "bottomleft", ext = c(-180, 180, -60, 90)) 
{
  length(dat)
  length(world)
  length(rst)
  fieldsmissing(dat, fields = c("ID", "x", "y", "x_original", 
                                "y_original", "Correction", "Exclude"))
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
  
  n <- length(x)
  sel <- rep(FALSE, length(x1))
  while(sum(sel) < n) {
    
    y <- coord2numeric(dat$y[ff])
    datid <- dat$ID[ff]
    xy1 <- data.frame(x = x, y = y)
    xy <- SpatialPoints(xy1)
    ex <- getextent(x, y, ext)
    xlm <- ex$xlm
    ylm <- ex$ylm
    plot(world, border = "gray", xlim = xlm, ylim = ylm, axes=T)
    dat$Correction <- as.character(dat$Correction)
    mtext(g1, side = 3, line = 3)
    vals <- extract(rst, xy)
    f2 <- which(!is.na(vals) & !sel)
    f3 <- which(is.na(vals) & !sel)
    if (nchar(group2) == 0) {
      points(x[f2], y[f2], col = "blue")
    }
    points(x[f3], y[f3], pch = 20, cex = 0.9, col = "red")
    if (nchar(group2) > 0) {
      f4 <- match(group2, cn)
      g2 <- unique(dat[, f4])
      ng2 <- length(g2)
      if (ng2 > 0) {
        cols <- rep(1, ng2)
        for (j in 1:ng2) {
          fg2 <- which(dat[ff, f4] == g2[j])
          points(x[fg2], y[fg2], col = j + 2)
          cols[j] <- j + 2
        }
        legend("bottomleft", legend = as.character(g2), col = cols, 
               pch = 18)
      }
    }
    fs <- identify(x[!sel], y[!sel], n = 1, plot = FALSE)
    if(!length(fs)) break
    fs <- which(!sel)[fs]
    idf <- datid[fs]
    ffa <- which(dat$ID == idf)
    idp <- dat$ID[ffa]
    if (locality == "") {
    }
    else {
      f1 <- match(locality, cn)
      mtext(dat[ffa, f1], side = 3, line = 0)
    }
    x2 <- x[fs]
    y2 <- y[fs]
    text(x2, y2, labels = idp, pos = 4, adj = 1)
    mtext(paste("x = ", x2, "; y = ", y2, sep = ""), side = 3, 
          line = 1, cex = 0.9)
    points(y2, x2, pch = 20, col = "purple")
    segments(x2, y2, y2, x2, lty = 2)
    points(x2 * -1, y2, pch = 20, col = "purple")
    segments(x2, y2, x2 * -1, y2, lty = 2)
    points(x2, y2 * -1, pch = 20, col = "purple")
    segments(x2, y2, x2, y2 * -1, lty = 2)
    points(x2 * -1, y2 * -1, pch = 20, col = "purple")
    segments(x2, y2, x2 * -1, y2 * -1, lty = 2)
    xs <- substddmm(x2)
    ys <- substddmm(y2)
    points(xs, ys, pch = 20, col = "purple")
    segments(x2, y2, xs, ys, lty = 2)
    signx <- ifelse(x2 < 0, -1, 1)
    signy <- ifelse(y2 < 0, -1, 1)
    xs2 <- abs(y2) * signx
    ys2 <- abs(x2) * signy
    points(xs2, ys2, pch = 20, col = "purple")
    segments(x2, y2, xs2, ys2, lty = 2)
    xnew <- c(y2, x2 * -1, x2, x2 * -1, xs, xs2, x2)
    ynew <- c(x2, y2, y2 * -1, y2 * -1, ys, ys2, y2)
    f1 <- identify(xnew, ynew, n = 1, plot = F)
    if (length(f1) == 0) 
      stop("Stopped")
    x1 <- xnew[f1]
    y1 <- ynew[f1]
    points(x1, y1, col = "black")
    mod <- format(Sys.time(), "%d-%m-%Y %H:%M:%S")
    dat$Modified <- as.character(dat$Modified)
    if (f1 == 7) {
      dat$Exclude[ffa] <- 1
      dat$Modified[ffa] <- mod
    }
    else {
      dc <- dat$Correction[ffa]
      dc <- str_replace(dc, "[.]", as.character(f1))
      dat$Correction[ffa] <- dc
      dat$x[ffa] <- x1
      dat$y[ffa] <- y1
      dat$x_original[ffa] <- x2
      dat$y_original[ffa] <- y2
      dat$Modified[ffa] <- mod
      #xyp <- paste("new x = ", dat$x[ff], "; y = ", dat$y[ff], 
      #             sep = "")
      #mtext(xyp, side = 1, line = 0, cex = 0.8)
      xyp <- paste("new x = ", dat$x[ffa], "; y = ", dat$y[ffa], "; Locality = ", dat$LocalityName[ffa], sep = "")
      print(xyp)
    }
  sel[fs] <- TRUE
  }
  return(dat)
}

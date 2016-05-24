alternatives <-
function (dat, group1 = "Species", group2 = "", world, rst, locality = "", 
          pos = "bottomleft", ext = c(-180, 180, -60, 90)) 
{
  length(dat)
  length(world)
  length(rst)
  fieldsmissing(dat, fields = c("ID", "x", "y", "x_original", 
                                "y_original", "Correction", "Modified", "Exclude"))

  dat$Correction <- as.character(dat$Correction)
  rn <- 1:nrow(dat)
  dat <- data.frame(rn, dat, stringsAsFactors = F)
  dat$Modified <- as.character(dat$Modified)
  cn <- names(dat)
  f3 <- match(group1, cn)
  ff <- which(dat$Exclude == 0)
  dat01 <- dat[ff, ]
  x1 <- coord2numeric(dat01$x)    
  y1 <- coord2numeric(dat01$y)
  nas <- cbind(is.na(x1) * 1, is.na(y1) * 1)
  f <- which(nas >= 1)
  dat01$Exclude[f] <- 1
  ff <- which(dat01$Exclude == 0)
  dat02 <- dat01[ff, ]
  x <- coord2numeric(dat02$x)
  
  n <- length(x)
  sel <- rep(FALSE, length(x1))
  while(sum(sel) < n) {
  
    y <- coord2numeric(dat02$y)
    datid <- dat02$ID
    xy1 <- data.frame(x = x, y = y)
    xy <- SpatialPoints(xy1)
    ex <- getextent(x, y, ext)
    xlm <- ex$xlm
    ylm <- ex$ylm
    plot(world, border = "gray", xlim = xlm, ylim = ylm, axes=T)
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
    idp <- dat02$ID[fs]
    x2 <- dat02$x[fs]
    y2 <- dat02$y[fs]
    ffa <- which(dat02$x == x2 & dat02$y == y2)
    spe <- dat02$Species[fs]
    fs2 <- which(dat02$Species == spe)
    points(dat02$x[fs2], dat02$y[fs2], col = "black", pch = 18)
    mtext(spe, side = 3, line = 3)
    if (locality == "") {
    }
    else {
      f1 <- match(locality, cn)
      mtext(dat02[fs, f1], side = 3, line = 0)
    }
    text(x2, y2, labels = idp, pos = 4, adj = 1)
    mtext(paste("x = ", x2, "; y = ", y2, sep = ""), side = 3, 
          line = 1, cex = 0.9)
    xs <- substddmm(x2)
    ys <- substddmm(y2)
    signx <- ifelse(x2 < 0, -1, 1)
    signy <- ifelse(y2 < 0, -1, 1)
    xs2 <- abs(y2) * signx
    ys2 <- abs(x2) * signy
    xnew <- c(y2, x2 * -1, x2, x2 * -1, xs, xs2, x2)
    ynew <- c(x2, y2, y2 * -1, y2 * -1, ys, ys2, y2)
    ni <- length(xnew) - 1
    for (j in 1:ni) {
      points(xnew[j], ynew[j], pch = 20, col = "purple")
      segments(x2, y2, xnew[j], ynew[j], lty = 2)
    }
    f1 <- identify(xnew, ynew, n = 1, plot = F)
    mod <- format(Sys.time(), "%d-%m-%Y %H:%M:%S")
    dat02$Modified <- as.character(dat02$Modified)
    if (length(f1) == 0) 
      stop("Stopped")
    x1 <- xnew[f1]
    y1 <- ynew[f1]
    points(x1, y1, col = "black")
    if (f1 == 7) {
      dat02$Exclude[ffa] <- 1
      dat02$Modified[ffa] <- mod
    } else {
      dc <- dat02$Correction[f1]
      dat02$Modified[ffa] <- mod
      dc <- str_replace(dc, "[.]", as.character(f1))
      dat02$Correction[ffa] <- dc
      dat02$x[ffa] <- x1
      dat02$y[ffa] <- y1
      dat02$x_original[ffa] <- x2
      dat02$y_original[ffa] <- y2
      xyp <- paste("new x = ", x1, "; y = ", y1, "; Locality = ", dat02$LocalityName[ffa], sep = "")
      print(xyp)
    }
    fi <- match(dat02$rn, dat$rn)
    dat[fi, ] <- dat02
    sel[fs] <- TRUE
  }
  dat <- dat[, -1]
  return(dat)
}

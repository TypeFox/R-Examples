persp.pemt <-
function(x, main, mar, ask = TRUE, col = "white", ...) {

  # 3D-Plot transition probabilities matrices 2D
  # through no ellispoidal interpolation
  #
  #          x object return by the function pemt()
  #       main title string
  #        mar vector to pass to par()  
  #        ask boolean to pass to par()
  #        col vector of colors ordered from 1 to 0
  #        ... other args to pass to image()

  ix <- length(x)
  nk <- x[[ix]]$nk
  mpoints <- x[[ix]]$mpoints
  coordsnames <- x[[ix]]$coordsnames
  nomi <- x[[ix]]$nomi

  if (missing(main) || !is.character(main)) main <- "Multidirectional transiogram"
  which.dire <- x[[ix]]$which.dire
  nimg <- dim(which.dire)

  ly <- matrix(c(rep(1, nk), nk^2+(1:nk)+2, 2:(nk^2+1), rep(0, nk), rep(nk^2+2, nk)), ncol = nk, byrow = TRUE)
  yl <- c(rep(0, nk + 4))
  yl1 <- c(rep(0, 2), (nk+1)*nk+2+(1:nk), rep(0, 2))
  ly <- cbind(yl, yl1, ly, yl, yl)
  widths <- c(rep(4 - 0.5 * nk, 2), rep(30/nk, nk), rep(4 - 0.5 * nk, 2)) / 4
  heights <- c(0.75, 1/3.5, rep(7.5/nk, nk), 2/3.5, 1)
  ly <- layout(ly, widths = widths, heights = heights, respect = TRUE)
  
  if (missing(mar) || !is.numeric(mar)) {
    mar <- sum(c(0, 1, 0, -0.5, -0.2, rep(0, 4), 0.2, rep(c(-0.1, rep(0, 3)), 2), 0)[1:nk])
    mar <- rep(mar, 4)
  }
  else {
    if (length(mar) < 4) {
      mar <- rep(mar[1], 4)
    }
    else {
      mar <- mar[1:4]
    }
  }

  oldcontour <- contour
  storage.mode(which.dire) <- "integer"

  for (i in 1:nimg[2]) {
    if (ask && (nimg[2] > 1)) {
      devAskNewPage(TRUE)
      on.exit(devAskNewPage(FALSE))
    }

    par(mar = c(0.3, 0.1, 0.3, 0.1))
    plot.new()
    title <- paste(main, " (", coordsnames[which.dire[1, i]], ", ", coordsnames[which.dire[2, i]], ")", sep = "")
    text(0.5, 0.5, labels = title, cex = 2)

    par(mar = mar)
    for (j in 1:nk) {
      for (k in 1:nk) {
        zzz <- t(matrix(rev(x[[i]]$Eprobs[j, k, ]), mpoints, mpoints))
        colid <- 1
        if (length(col) > 1) {
          colth <- seq(1 / length(col), 1, length.out = length(col))
          colid <- sapply(zzz, function(xx) sum(xx < colth))
          dim(colid) <- dim(zzz)
          colid <- colid[-1, -1]
          zlim <- range(zzz, na.rm = TRUE)
        }
        else {
          zlim <- 0:1
        }
        persp(x[[i]]$X, x[[i]]$Y, zzz, 
              col = rev(col)[colid], zlab = "", zlim = zlim,
              xlab = coordsnames[which.dire[1, i]], 
              ylab = coordsnames[which.dire[2, i]], ...) #sub = paste(nomi[j], "to", nomi[k])
      }
    }

    par(mar = c(2.2, 0.1, 1.2, 0.1))
    if (length(col) > 1) {
      image(matrix(colth, ncol = 1), col = col, main="Legend", axes = FALSE)
      box()
      axis(1)
    }
    else {
      plot.new()
    }
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    for (j in c(0, 90)) {
      for (lbls in nomi) {
        plot.new()
        text(0.5, 0.5, labels = lbls, srt = j, font = 3)
      }
    }
  }

}

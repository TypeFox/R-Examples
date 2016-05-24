image.pemt <-
function(x, main, mar, ask = TRUE, ..., nlevels = 10, contour = FALSE) {
  # Plot transition probabilities matrices 2D
  # through no ellispoidal interpolation
  #
  #          x object return by the function pemt()
  #       main title string
  #        mar vector to pass to par()  
  #        ask boolean to pass to par()
  #        ... other args to pass to image()
  #    nlevels number of levels to pass to contour()
  #    contour boolean values that permits to draw contour lines

  ix <- length(x)
  nk <- x[[ix]]$nk
  mpoints <- x[[ix]]$mpoints
  coordsnames <- x[[ix]]$coordsnames

  if (missing(main) || !is.character(main)) main <- "Pseudoempirical transiogram"
  which.dire <- x[[ix]]$which.dire
  nimg <- dim(which.dire)

  ly <- matrix(c(rep(1, nk), rep(0, nk), 2:(nk^2+1), rep(0, nk), rep(nk^2+2, nk)), ncol = nk, byrow = TRUE)
  yl <- c(rep(0, nk + 4))
  ly <- cbind(yl, yl, ly, yl, yl)
  widths <- c(rep(4 - 0.5 * nk, 2), rep(30/nk, nk), rep(4 - 0.5 * nk, 2)) / 4
  heights <- c(0.75, 1/3.5, rep(7.5/nk, nk), 2/3.5, 1)
  ly <- layout(ly, widths = widths, heights = heights, respect = TRUE)
  
  nomi <- x[[ix]]$nomi
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
        image(x[[i]]$X, x[[i]]$Y,
              t(matrix(rev(x[[i]]$Eprobs[j, k, ]), mpoints, mpoints)),
              xlab = paste("Coord", which.dire[1, i]),
              ylab = paste("Coord", which.dire[2, i]), ..., 
              zlim = 0:1, axes = FALSE)
        box()
        if (k == 1) axis(2)
        if (k == nk) axis(4, 0, labels = nomi[j], tick = FALSE, font = 3)
        if (j == nk) axis(1)
        if (j == 1) axis(3, 0, labels = nomi[k], tick = FALSE, font = 3)
        if (is.null(x[[i]]$Tprobs)) {
          message("Contour values missing or not computable")
        }
        else {
          if (contour) {
            contour(x[[i]]$X, x[[i]]$Y,
                    t(matrix(rev(x[[i]]$Tprobs$Tmat[j, k, ]), mpoints, mpoints)),
                    nlevels, add = TRUE)
          }
        }
      }
    }

    par(mar = c(2.2, 0.1, 1.2, 0.1))
    image(matrix(0:500/500, ncol=1), ..., main="Legend", axes = FALSE)
    box()
    axis(1)

  }

}

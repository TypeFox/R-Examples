niche.plot <-
function(niche.par, niche.data, alpha = .95,
                       species.names, iso.names,
                       col, ndens = 512, pfrac = 0, xlab) {
  niso <- ncol(niche.par[[1]]$mu)
  nspec <- length(niche.par)
  npts <- 100 # number of points for each ellipse
  nell <- sapply(niche.par, function(x) nrow(x$mu)) # number of ellipses per species
  if(missing(species.names)) species.names <- names(niche.par)
  if(missing(iso.names)) iso.names <- colnames(niche.par[[1]]$mu)
  # create all the ellipses to get the plot limits right.
  ell <- vector("list", nspec)
  names(ell) <- names(species.names)
  D <- combn(niso, 2)
  for(ii in 1:nspec) {
    ell.tmp <- array(NA, c(nell[ii], ncol(D), npts+1, 2))
    for(jj in 1:nell[ii]) {
      for(kk in 1:ncol(D)) {
        ell.coord <- ellipse(niche.par[[ii]]$mu[jj, D[,kk]],
                             V = niche.par[[ii]]$Sigma[D[,kk], D[,kk], jj],
                             alpha = alpha, n = npts)
        ell.tmp[jj,kk,,] <- ell.coord
      }
    }
    ell[[ii]] <- ell.tmp
  }
  # plot limits.
  lims <- array(sapply(niche.data, function(x) apply(x, 2, range)),
                dim = c(2, niso, nspec))
  lims <- apply(lims, 2, range)
  # plots
  par(mfcol = c(niso,niso), mar = rep(.5, 4), oma = rep(4,4))
  for(ci in 1:niso) {
    for(ri in 1:niso) {
      # initialize plot
      plot.new()
      plot.window(lims[,ci], lims[,ri])
      if (ci == ri) {
        # diagonals: density plots
        xdens <- matrix(NA, ndens, nspec)
        ydens <- xdens
        for(ii in 1:nspec) {
          den <- density(niche.data[[ii]][,ci], n = ndens)
          xdens[,ii] <- den$x
          ydens[,ii] <- den$y
        }
        for(ii in 1:nspec) {
          ly <- par("usr")[1:2]
          ly[2] <- ly[1] + pfrac*(ly[2]-ly[1])
          ly[3] <- (ly[2]-ly[1])/nspec
          segments(x0 = niche.data[[ii]][,ci],
                   y0 = ly[1]+(ii-1)*ly[3], y1 = ly[1]+ii*ly[3], col = col[ii])
          ly <- ly[2] + ydens[,ii]/max(ydens)*(lims[2,ci]-ly[2])
          lines(xdens[,ii], ly, col = col[ii])
        }
      }
      if (ri > ci) {
        # lower triangle: point plots
        for(ii in 1:nspec) {
          points(niche.data[[ii]][,c(ci,ri)], col = col[ii], pch = 16)
        }
      }
      if (ri < ci) {
        # upper triangle: ellipses
        for(ii in 1:nspec) {
          for(jj in 1:nell[ii]) {
            lines(ell[[ii]][jj,which(D[1,] == ri & D[2,] == ci),,2:1], col = col[ii])
          }
        }
      }
      box()
      if(ci == niso) axis(side = 4) else axis(side = 4, labels = FALSE)
      if(ri == niso) axis(side = 1) else axis(side = 1, labels = FALSE)
      if(ri == 1) mtext(text = iso.names[ci], side = 3, line = 1)
      if(ci == 1) mtext(text = iso.names[ri], side = 2, line = 1)
    }
  }
  if(!missing(xlab)) {
    mtext(text = xlab, side = 1, outer = TRUE, line = 2.2, cex = .9)
  }
  legend(x = "topleft", legend = species.names, fill = col, bty = "n", cex = 1.25)
}

PlotEquiMap <- function(var, lon, lat, toptitle = '', sizetit = 1, units = '', 
                        brks = NULL, cols = NULL, square = TRUE, 
                        filled.continents = TRUE, contours = NULL, 
                        brks2 = NULL, dots = NULL, axelab = TRUE, labW = FALSE, 
                        intylat = 20, intxlon = 20, drawleg = TRUE, 
                        boxlim = NULL, boxcol = 'purple2', boxlwd = 10,
                        subsampleg = 1, numbfig = 1, colNA = 'white') {
  #library(GEOmap)
  #library(geomapdata)
  #library(maps)
  data(coastmap, package = 'GEOmap', envir = environment())
  #
  #  Input arguments 
  # ~~~~~~~~~~~~~~~~~
  #
  dims <- dim(var)
  if (length(dims) > 2) {
    stop("Only 2 dimensions expected for var : (lon,lat) ")
  }
  if (dims[1] != length(lon) | dims[2] != length(lat)) {
    if (dims[1] == length(lat) & dims[2] == length(lon)) {
      var <- t(var)
      dims <- dim(var)
    } else {
      stop("Inconsistent var dimensions / longitudes +  latitudes")
    }
  }
  latb <- sort(lat, index.return = TRUE)
  dlon <- lon[2:dims[1]] - lon[1:(dims[1] - 1)]
  wher <- which(dlon > (mean(dlon) + 1))
  if (length(wher) > 0) {
    lon[(wher + 1):dims[1]] <- lon[(wher + 1):dims[1]] - 360
  }
  lonb <- sort(lon, index.return = TRUE)
  latmin <- floor(min(lat) / 10) * 10
  latmax <- ceiling(max(lat) / 10) * 10
  lonmin <- floor(min(lon) / 10) * 10
  lonmax <- ceiling(max(lon) / 10) * 10

  colorbar <- colorRampPalette(c("dodgerblue4", "dodgerblue1",
              "forestgreen", "yellowgreen", "white", "white",
              "yellow", "orange", "red", "saddlebrown"))
  if (is.null(brks) == TRUE) {
    ll <- signif(min(var, na.rm = TRUE), 4)
    ul <- signif(max(var, na.rm = TRUE), 4)
    if (is.null(cols) == TRUE) {
      cols <- colorbar(10)
    }
    nlev <- length(cols)
    brks <- signif(seq(ll, ul, length.out = 1 + nlev), 4)
  } else {
    if (is.null(cols) == TRUE) {
      nlev <- length(brks) - 1
      cols <- colorbar(nlev)
    } else {
      if (length(cols) != (length(brks) - 1)) {
        stop("Inconsistent colour levels / list of colours")
      }
    }
  }

  if (is.null(brks2) == TRUE ) {
    if (is.null(contours)) { 
      if (square == FALSE) {
        brks2 <- brks
        contours <- var 
      }
    } else {
      ll <- signif(min(contours, na.rm = TRUE), 2)
      ul <- signif(max(contours, na.rm = TRUE), 2)
      brks2 <- signif(seq(ll, ul, length.out = length(brks)), 2)
    }
  }
  #
  #  Plotting the map
  # ~~~~~~~~~~~~~~~~~~
  #
  if (axelab == TRUE) { 
    ypos <- seq(latmin, latmax, intylat)
    xpos <- seq(lonmin, lonmax, intxlon)
    letters <- array('', length(ypos))
    letters[ypos < 0] <- 'S'
    letters[ypos > 0] <- 'N'
    ylabs <- paste(as.character(abs(ypos)), letters, sep = '')
    letters <- array('', length(xpos))
    if (labW) {
      nlon <- length(xpos)
      xpos2 <- xpos  
      xpos2[xpos2 > 180] <- 360 - xpos2[xpos2 > 180]
    }
    letters[xpos < 0] <- 'W'
    letters[xpos > 0] <- 'E'
    if (labW) {
      letters[xpos == 0] <- ' '
      letters[xpos == 180] <- ' '
      letters[xpos > 180] <- 'W'  
      xlabs <- paste(as.character(abs(xpos2)), letters, sep = '')
    } 
    else{
      xlabs <- paste(as.character(abs(xpos)), letters, sep = '')
    }
    xmargin <- 1.2 - (numbfig ^ 0.2) * 0.05
    ymargin <- 3 - (numbfig ^ 0.2)
    spaceticklab <- 1.3 - (numbfig ^ 0.2) * 0.8
    topmargin <- 0.4
    ymargin2 <- 1.5 - (numbfig ^ 0.2) * 0.9
  } else {
    xmargin <- 0.2
    ymargin <- 0.2 + switch(as.character(square), 'FALSE' = 1.8, 0)
    topmargin <- 0.2
    spaceticklab <- 1
    ymargin2 <- 0.2
  }
  if (toptitle != '') {
    topmargin <- 2.5 - (numbfig ^ 0.2) * 0.6
  }
  if (min(lon) < 0) {
    continents <- 'world'
  } else {
    continents <- 'world2'
  }
  if (square) {
    if (drawleg) {
      layout(matrix(1:2, ncol = 1, nrow = 2), heights = c(5, 1))
    }
    par(mar = c(xmargin, ymargin, topmargin, ymargin2), cex = 1.4, 
        mgp = c(3, spaceticklab, 0), las = 0)
    if (colNA != 'white') {
      blanks <- array(0, dim = c(length(lonb$x), length(latb$x)))
      image(lonb$x, latb$x, blanks, col = c(colNA), breaks = c(-1, 1), 
            main = toptitle, cex.main = (1.5 / numbfig ^ (0.2)) * sizetit, 
            axes = FALSE, xlab = "", ylab = "")
      flagadd <- T
    } else {  
      flagadd <- F
    }
    image(lonb$x, latb$x, var[lonb$ix, latb$ix], col = cols, breaks = brks, 
          main = toptitle, axes = FALSE, xlab = "", ylab = "", 
          cex.main = (1.5 / numbfig ^ (0.2)) * sizetit, add = flagadd)
    if (axelab == TRUE) {
      axis(2, at = ypos, labels = ylabs, cex.axis = 1 / (numbfig ^ 0.3), 
           tck = -0.03)
      axis(1, at = xpos, labels = xlabs, cex.axis = 1 / (numbfig ^ 0.3),
           tck = -0.03)
    }
    if (is.null(contours) == FALSE) {
      contour(lonb$x, latb$x, contours[lonb$ix, latb$ix], levels = brks2,
              method = "edge", add = TRUE, labcex = 1 / numbfig, 
              lwd = 0.5 / (numbfig ^ 0.5))
    }
    map(continents, interior = FALSE, add = TRUE, lwd = 1)
    box()
  } else {
    par(mar = c(xmargin + 5, ymargin + 1.5, topmargin, ymargin2), 
        cex.main = (1.6 * numbfig ^ (0.3)) * sizetit, cex.axis = 1.4, 
        cex.lab = 1.6, mgp = c(3, spaceticklab + 0.5, 0), las = 0)
    if (axelab == TRUE) {
      filled.contour(lonb$x, latb$x, var[lonb$ix, latb$ix], xlab = "", 
                     levels = brks, col = cols, ylab = "", main = toptitle,
                     key.axes = axis(4, brks[seq(1, length(brks), subsampleg)],
                     cex.axis = 1 / (numbfig ^ 0.3)), 
                     plot.axes = {axis(2, at = ypos, labels = ylabs, 
                                  cex.axis = 1 / (numbfig ^ 0.3), tck = -0.03)
                                  axis(1, at = xpos, labels = xlabs,
                                  cex.axis = 1 / (numbfig ^ 0.3), tck = -0.03)
                                  contour(lonb$x, latb$x, contours[lonb$ix,
                                          latb$ix], levels = brks2, 
                                          method = "edge", add = TRUE,
                                          labcex = 1, lwd = 2)
                                  map(continents, interior = FALSE, 
                                      xlim = c(lonmin, lonmax),
                                      ylim = c(latmin, latmax), add = TRUE)},
                     key.title = title(main = units, 
                                 cex.main = (1.2 * numbfig ^ (0.3)) * sizetit))
    } else {
      filled.contour(lonb$x, latb$x, var[lonb$ix, latb$ix], xlab = "",
                     levels = brks, col = cols, ylab = "", main = toptitle,
                     key.axes = axis(4, brks[seq(1, length(brks), subsampleg)],
                     cex.axis = 1 / (numbfig ^ 0.3)), 
                     plot.axes = {contour(lonb$x, latb$x, contours[lonb$ix, 
                                          latb$ix], levels = brks2, 
                                          method = "edge", add = TRUE,
                                          labcex = 1, lwd = 2)
                                  map(continents, interior = FALSE, 
                                      xlim = c(lonmin, lonmax), 
                                      ylim = c(latmin, latmax), add = TRUE)},
                     key.title = title(main = units, 
                                 cex.main = (1.2 * numbfig ^ (0.3)) * sizetit))
    }
  } 
  #
  #  Adding black dots  
  # ~~~~~~~~~~~~~~~~~~~
  #
  if (is.null(dots) == FALSE) {
    for (ix in 1:length(lon)) {
      for (jy in 1:length(lat)) {
        if (is.na(var[ix, jy]) == FALSE) {
          if (dots[ix, jy] == TRUE) {
            text(x = lon[ix], y = lat[jy], ".", 
                 cex = 12 / (sqrt(sqrt(length(var))) * numbfig ^ 0.5))
          }
        }
      }
    }
  }
  #
  #  Filling continents
  # ~~~~~~~~~~~~~~~~~~~~
  # 
  if (square == TRUE & filled.continents == TRUE) {
    if (min(lon) >= 0) {
      ylat <- latmin:latmax
      xlon <- lonmin:lonmax
      proj <- setPROJ(1, LON0 = mean(xlon), LAT0 = mean(ylat), LATS = ylat,
                      LONS = xlon)
      coastmap$STROKES$col[which(coastmap$STROKES$col == "blue")] <- "white"
      par(new = TRUE)
      plotGEOmap(coastmap, PROJ = proj, border = 'black', add = TRUE)
      box()
    } else {
      map(continents, interior = FALSE, wrap = TRUE, lwd = 0.7, col = gray(0.5),
          fill = TRUE, add = TRUE, border = gray(0.5))
    }
  }
  # Draw rectangle on the map
  if (is.null(boxlim) == FALSE) {
    boxlimaux <- boxlim
    if(boxlim[1]>boxlim[3]){
      boxlimaux[1] <- boxlim[1]-360
    }
    if (length(boxlimaux)!=4){
      stop('Region to be highlighted is ill defined')
    } else if(boxlimaux[2] < latmin | boxlimaux[4] > latmax | boxlimaux[1]<lonmin | boxlimaux[3]>lonmax){
      stop('Invalid boundaries')
    } else if(boxlimaux[1]<0 && boxlimaux[3]>0){
      #segments south
      segments(boxlimaux[1], boxlimaux[2],0, boxlimaux[2], col=boxcol, lwd=boxlwd)
      segments(0,boxlimaux[2], boxlimaux[3], boxlimaux[2], col=boxcol, lwd=boxlwd) 
      #segments north
      segments(boxlimaux[1], boxlimaux[4],0, boxlimaux[4], col=boxcol, lwd=boxlwd)
      segments(0,boxlimaux[4], boxlimaux[3], boxlimaux[4], col=boxcol, lwd=boxlwd) 
      #segments west
      segments(boxlimaux[1], boxlimaux[2], boxlimaux[1],boxlimaux[4], col=boxcol, lwd=boxlwd )  
      #segments est
      segments(boxlimaux[3], boxlimaux[2], boxlimaux[3],boxlimaux[4], col=boxcol, lwd=boxlwd )          
    } else {
      rect(boxlimaux[1], boxlimaux[2], boxlimaux[3], boxlimaux[4], border=boxcol, col=NULL, lwd=boxlwd, lty='solid')              
    }
  }
  #
  #  Colorbar
  # ~~~~~~~~~~
  #
  if (square & drawleg) {
    par(mar = c(1.5, ymargin + 1.5, 2.5, ymargin2), mgp = c(1.5, 0.3, 0), 
        las = 1, cex = 1.2)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols, 
          xlab = '', ylab = '', main = units, cex.main = 1.1)
    box()
    axis(1, at = seq(0.5, length(brks) - 0.5, subsampleg), 
         labels = brks[seq(1, length(brks), subsampleg)])
    #title(xlab = units)
  }
}

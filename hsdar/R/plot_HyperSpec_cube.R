.plot_side_cube <- function(ra, fi, side, col = colorRamp(palette(heat.colors(100))), ncol = ncol, nrow = nrow)
{ 
  if (side == 1)
  {
    test <- spectra(getValuesBlock(ra, row = nrow[1], nrows = 1, col = 1))
    png_height <- nrow(ra)
  }
  if (side == 2)
  {
    test <- spectra(getValuesBlock(ra, row = 1, nrows = nrow(ra), col = ncol[1], ncol = 1))
    png_height <- ncol(ra)
  }
  if (side == 3)
  {
    test <- spectra(getValuesBlock(ra, row = 1 + nrow(ra) - nrow[2], nrows = 1, col = 1))
    png_height <- nrow(ra)
  }
  if (side == 4)
  {
    test <- spectra(getValuesBlock(ra, row = 1, nrows = nrow(ra), col = 1 + ncol(ra) - ncol[2], ncol = 1))
    png_height <- ncol(ra)
  }
  
  png(filename = fi, height = png_height, width = nbands(ra))
 
  test <- (test - min(test))/(max(test)-min(test))
  test_col <- rgb(col(test), maxColorValue = 255)
  test_col <- matrix(test_col, ncol = ncol(test), nrow = nrow(test))
  
  par(mar = c(0,0,0,0), oma = c(0,0,0,0))
  plot(c(0,1), c(0,1), type = "l", yaxt = "n", bty = "n", xaxs = "i", yaxs = "i")
  rasterImage(test_col, 0, 0, 1, 1, interpolate = TRUE)
  dev.off()
}

.plotRGB_temp <- function(ras, side = NULL, ncol = 1, nrow = 1, sidecol = colorRamp(palette(heat.colors(100))), ...)
{
  fi <- tempfile(fileext = ".png")
  if (is.null(side))
  {
    png(filename = fi, width = ncol(ras), height = nrow(ras))
    plotRGB(ras, ...)
    dev.off()
  } else {
    if (length(ncol) == 1)
      ncol <- c(ncol, ncol)
    if (length(nrow) == 1)
      nrow <- c(nrow, nrow)  
      
    .plot_side_cube(ra = ras, fi = fi, side = side, ncol = ncol, nrow = nrow, col = sidecol)
  }    
  return(fi)
}

cubePlot <- function(x, r, g, b, ncol = 1, nrow = 1, sidecol = colorRamp(palette(heat.colors(100))), ...)
{
  if (!requireNamespace("rgl", quietly = TRUE))
    stop("Library 'rgl' is required to plot 3D cube of hyperspectral data")
  
  ra <- x
  
  if (missing(r))
    r <- which.min(abs(wavelength(ra) - 680))
    
  if (missing(g))
    g <- which.min(abs(wavelength(ra) - 540))
    
  if (missing(b))
    b <- which.min(abs(wavelength(ra) - 470))
  
  z <- matrix(0, ncol = ncol(ra), nrow = nrow(ra))
  x <- c(1:nrow(z))-1
  y <- c(1:ncol(z))-1

  texture <- .plotRGB_temp(ras = ra, r = r, g = g, b = b, stretch = "hist")
  texture_side_1 <- .plotRGB_temp(ras = ra, side = 1, ncol = ncol, nrow = nrow, sidecol = sidecol)
  texture_side_2 <- .plotRGB_temp(ras = ra, side = 2, ncol = ncol, nrow = nrow, sidecol = sidecol)
  texture_side_3 <- .plotRGB_temp(ras = ra, side = 3, ncol = ncol, nrow = nrow, sidecol = sidecol)
  texture_side_4 <- .plotRGB_temp(ras = ra, side = 4, ncol = ncol, nrow = nrow, sidecol = sidecol)

  rgl::open3d()
  rgl::rgl.surface(x, y, z, col = "white", 
                   axes = FALSE, box = FALSE,
                   texture = texture)

  z <- matrix(0, ncol = nrow(ra), nrow = nbands(ra))
  x <- c(1:nrow(z))-1
  y <- c(1:ncol(z))-1
  rgl::rgl.surface(x, y, z, col = "white",
                   axes = FALSE, box = FALSE, coords = c(2,3,1),
                   texture = texture_side_1, add = TRUE)

  z <- z + ncol(ra)-1
  rgl::rgl.surface(x, y, z, col = "white",
                   axes = FALSE, box = FALSE, coords = c(2,3,1),
                   texture = texture_side_3, add = TRUE)

  z <- matrix(0, ncol = ncol(ra), nrow = nbands(ra))
  x <- c(1:nrow(z))-1
  y <- c(1:ncol(z))-1
  rgl::rgl.surface(x, y, z, col = "white",
                   axes = FALSE, box = FALSE, coords = c(2,1,3),
                   texture = texture_side_2, add= TRUE)

  z <- z + nrow(ra)-1
  rgl::rgl.surface(x, y, z, col = "white",
                   axes = FALSE, box = FALSE, coords = c(2,1,3),
                   texture = texture_side_4, add= TRUE)
}


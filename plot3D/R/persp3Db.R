## =============================================================================
## main function, with input of colors, no colorvar
## =============================================================================

persp3Db <- function(x = seq(0, 1, length.out = nrow(col) +1), 
                     y = seq(0, 1, length.out = ncol(col) +1), 
                     z, col, ..., 
                     phi = 40, theta = 40, NAcol = "white", breaks = NULL,
                     border = NA, facets = TRUE, panel.first = NULL, bty = "b", 
                     lighting = FALSE, shade = NA, ltheta = -135, lphi = 0,
                     add = FALSE, plot = TRUE){


  if (is.vector(x)) {
    if (length(x) == nrow(col))
      x <- extendvec(x)
    if (length(x) != nrow(col) + 1)
      stop ("length of 'x' should be = rows of 'col' +1")

  } else  if (is.matrix(x)) {
    if (nrow(x) == nrow(col) & ncol(x) == ncol(col))
      x <- extend(x)
    if (nrow(x) != nrow(col) + 1)
      stop ("rows of 'x' should be = rows of 'col' +1")
    if (ncol(x) != ncol(col) + 1 )
      stop ("columns of 'x' should be = columns of 'col' + 1")
  }
    
  if (is.vector(y)) { 
    if (length(y) == ncol(col))
      y <- extendvec(y)
    if (length(y) != ncol(col) + 1)
      stop ("length of 'y' should be = columns of 'col' +1")
  } else if (is.matrix(y)) {
    if (nrow(y) == nrow(col) & ncol(y) == ncol(col))
      y <- extend(y)
    if (nrow(y) != nrow(col) + 1)
      stop ("rows of 'y' should be = rows of 'col' +1")
    if (ncol(y) != ncol(col) + 1 )
      stop ("columns of 'y' should be = columns of 'col' + 1")
  }
  
  if (nrow(z) == nrow(col) & ncol(z) == ncol(col))
    z <- extend(z)

  if (nrow(z) != nrow(col) + 1)
    stop ("rows of 'z' should be = rows of 'col' +1")
    
  if (ncol(z) != ncol(col) + 1 )
    stop ("columns of 'z' should be = columns of 'col' + 1")


  plist <- initplist(add)
        
  dot <- splitdotpersp(list(...), bty, lighting, 
    x, y, z, plist = plist, shade, lphi, ltheta, breaks = breaks)

  if (is.null(plist)) {
    do.call("perspbox", c(alist(x, y, z,  
                     phi = phi, theta = theta, plot = plot, 
                     colkey = FALSE, col = col), dot$persp))
    plist <- getplist()
  }
  
  if (is.function(panel.first)) 
    panel.first(plist$mat)         

  if (! is.matrix(x)) { 
    x <- matrix(nrow = nrow(z), ncol = ncol(z), data = x)
    y <- matrix(nrow = nrow(z), ncol = ncol(z), data = y, byrow = TRUE)
  } 
   
  lwd <- ifelse (is.null (dot$points$lwd), 1, dot$points$lwd)
  lty <- ifelse (is.null (dot$points$lty), 1, dot$points$lty)

  sl <- Sortlist(x, y, z, plist, Polar = FALSE)

  if (dot$shade$type != "none") 
    col <- facetcols (x, y, z, col, dot$shade, Extend = FALSE)
  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA

 # Draw colored polygons           
  Poly <- list()
  Poly$img <- list(list(x = x, y = y, z = z, col = col, sl = sl, 
    NAcol = NAcol, facets = facets, border = border, lwd = lwd, 
    lty = lty, alpha = alpha, mapped = FALSE))

 # plot it
  plist <- plot.struct.3D(plist, poly = Poly, plot = plot)  

  setplist(plist)   
  invisible(plist$mat)
}                  


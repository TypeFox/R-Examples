## =============================================================================
## 3-d representation using slices in x, y or z
## =============================================================================
# x, y, z vectors or arrays, colvar: array

addslice <- function(poly, x, y, z, colvar, xs = NULL, 
                     ys = NULL, zs = NULL, plist,
                     col = NULL, NAcol = "white", breaks = NULL,
                     border = NA, facets = TRUE, lwd = 1, lty = 1, 
                     clim = NULL, shadedot = NULL,
                     lighting = FALSE, alpha = NULL)  {

  if (! ispresent(colvar))
    stop("'colvar' has to be defined and be an array of dimension 3")

 # check dimensionality 
  DD <- dim(colvar)
  if (length(DD) != 3)
    stop("'colvar' has to be an array of dimension 3")
  if (DD[1] !=  length(x))
    stop("dimension of 'colvar' not compatible with length of 'x'")
  if (DD[2] !=  length(y))
    stop("dimension of 'colvar' not compatible with length of 'y'")
  if (DD[3] !=  length(z))
    stop("dimension of 'colvar' not compatible with length of 'z'")

  crange <- diff(clim)
  N      <- length(col) -1

 # Colors for NA values   
  if (any (is.na(colvar)) & ! is.null(NAcol) ) {
    CC <- checkcolors(colvar, col, NAcol, clim)
    clim   <- CC$lim
    col    <- CC$col
    colvar <- CC$colvar
  }
 
 # Function to add images on a plane as polygons
  image.plane <- function(xs, ys, zs, paint = FALSE, i = 0) {
    if (i == 1)   # all xs are same
      ix <-  FindInterval(xs[1], x, all.inside = TRUE)
    else      
      ix <- FindInterval(xs, x, all.inside = TRUE)
    if (i == 2) 
      iy <- FindInterval(ys[1], y, all.inside = TRUE)
    else 
      iy <- FindInterval(ys, y, all.inside = TRUE)
    if (i == 3) 
      iz <- FindInterval(zs[1], z, all.inside = TRUE)
    else 
      iz <- FindInterval(zs, z, all.inside = TRUE)
    
   # colorvar 
    cv <- matrix(nrow = nrow(xs), ncol = ncol(xs), data = colvar[cbind(ix, iy, iz)])

   # add polygon
    poly <<- addimg(poly, xs, ys, zs, colvar = cv, plist = plist, 
        col = col, NAcol = NAcol, breaks = breaks,
        border = border, facets = facets,
        resfac = 1, clim = clim, lwd = lwd, lty = lty, 
        ltheta = shadedot$ltheta, lphi = shadedot$lphi, 
        shade = shadedot$shade, lighting = lighting, alpha = alpha)
    
  } # end function imageplane

 # Function to first create a plane and then draw an image on it
  add.plane <- function(xs, ys, zs, i = 0) {

    M <- mesh(xs, ys, zs)
    image.plane (M$x[,,], M$y[,,], M$z[,,], i = i) # [,,] to make sure it is an array

  }  # end function addplane


  if (any(diff(c( is.matrix(xs), is.matrix(ys), is.matrix(zs)))) != 0)
    stop ("'xs', 'ys' and 'zs' should be a matrix if one of them is")

  if (is.matrix(xs)) {  # xs,.. are matrices defining the plane on which to plot 
    if (any (dim(xs) != dim(ys)))
      stop("'xs' and 'ys' should have same dimension")
    
    if (any (dim(xs) != dim(zs)))
      stop("'xs' and 'zs' should have same dimension")
    
    image.plane(xs, ys, zs, paint = TRUE)  
  
  } else { # xs, ys, zs define the positions in x, y, z on which to plot
    if (! is.null(xs))
      for (x.s in xs[!is.na(xs)])
        add.plane(x.s, y, z, 1)
     
    if (! is.null(ys))
      for (y.s in ys[!is.na(ys)]) 
        add.plane(x, y.s, z, 2)
    
    if (! is.null(zs))
      for (z.s in zs[!is.na(zs)]) 
        add.plane(x, y, z.s, 3)
  }
  return(poly)
}                     

## =============================================================================
## main slice3D function
## =============================================================================

slice3D <- function(x, y, z, colvar, ..., 
                    phi = 40, theta = 40, 
                    xs = min(x),
                    ys = max(y),
                    zs = min(z),
                    col = NULL, NAcol = "white", breaks = NULL,
                    border = NA, facets = TRUE, 
                    colkey = NULL, panel.first = NULL,
                    clim = NULL, clab = NULL, bty = "b",
                    lighting = FALSE, shade = NA, ltheta = -135, lphi = 0, 
                    add = FALSE, plot = TRUE) {

  plist <- initplist(add)

  dot <- splitdotpersp(list(...), bty, lighting, x, y, z, plist = plist,
    shade, lphi, ltheta, breaks = breaks)

    if (is.null(col))
      if (is.null(breaks))
        col <- jet.col(100)
      else
        col <- jet.col(length(breaks)-1)

  iscolkey <- is.colkey(colkey, col)
  if (iscolkey) 
    colkey <- check.colkey(colkey)

  breaks <- check.breaks(breaks, col)
  if (! is.null(dot$alpha)) col <- setalpha(col, dot$alpha)

  if (is.null(clim))
    clim <- range(colvar, na.rm = TRUE)     

  if (dot$clog) {
    colvar <- log(colvar)
    clim <- log(clim)
  }
        
  if (is.null(plist)) {
    do.call("perspbox", c(alist(x = range(x), y = range(y), 
             z = range(z, na.rm = TRUE),
             phi = phi, theta = theta, plot = plot, 
             colkey = colkey, col = col), dot$persp))
    plist <- getplist()
  }  
  if (is.function(panel.first)) 
    panel.first(plist$mat)         

  lwd <- ifelse(is.null(dot$points$lwd), 1, dot$points$lwd)
  dot$points$lwd <- NULL

  lty <- ifelse(is.null(dot$points$lty), 1, dot$points$lty)
  dot$points$lty <- NULL

  Poly <- addslice(NULL, x, y, z, colvar, xs = xs, ys = ys, zs = zs, plist = plist,
                   col = col, NAcol = NAcol, breaks = breaks,
                   border = border, facets = facets,
                   clim = clim, shadedot = dot$shade, lwd = lwd, lty = lty,
                   lighting = lighting, alpha = dot$alpha)
   
  if (iscolkey)  
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
      dot$clog, type = "slice3D", breaks = breaks)
 
  plist <- plot.struct.3D(plist, poly = Poly, plot = plot)  

  setplist(plist)
  invisible(plist$mat)
}
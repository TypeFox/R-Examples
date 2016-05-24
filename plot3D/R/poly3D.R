
# a function to create polygons from vectors

add.poly <- function(poly, x, y, z, cv, col, NAcol, breaks,
  clim, facets, border, lwd = 1, lty = 1, ...)  {

  nr <- nrow(x) - 1
  nc <- ncol(x) - 1
    
 # create polygons
  ix <- rep(1:nr, nc)
  iy <- as.vector(matrix(nrow = nr, ncol = nc, 
                     data = 1:nc, byrow =TRUE))
  xx <- x
  yy <- y
  zz <- z
    
 # the polygon coordinates
  PolyX <- rbind(xx[cbind(ix,     iy    )],
                 xx[cbind(ix + 1, iy    )],
                 xx[cbind(ix + 1, iy + 1)],
                 xx[cbind(ix,     iy + 1)], NA)
  PolyY <- rbind(yy[cbind(ix,     iy    )],
                 yy[cbind(ix + 1, iy    )],
                 yy[cbind(ix + 1, iy + 1)],
                 yy[cbind(ix,     iy + 1)], NA)
  PolyZ <- rbind(zz[cbind(ix,     iy    )],
                 zz[cbind(ix + 1, iy    )],
                 zz[cbind(ix + 1, iy + 1)],
                 zz[cbind(ix,     iy + 1)], NA)

 # colvar is converted to colors.
  if (! is.null(cv)) {
    if (length(cv) == nrow(x))    # take colvar averages
      cv <- 0.5*(cv[-1] + cv[-length(cv)])

 # Colors for values = NA
    if (any (is.na(cv)) & ! is.null(NAcol) ) {
      CC <- checkcolors(cv, col, NAcol, clim)
      clim   <- CC$lim
      col    <- CC$col
      cv     <- CC$colvar
    }
    if (is.null(breaks)) {
      crange <- diff(clim)
      N      <- length(col) -1
      Col <- col[1 + trunc((cv - clim[1])/crange*N)]
    } else {
      zi <- .bincode(cv, breaks, TRUE, TRUE)
      Col <- col[zi]
      Col[is.na(Col)] <- NAcol
    }
  } else 
    Col <- rep(col, length.out = length(cv))
      
  Lty <- rep(lty, length.out = length(cv))
  Lwd <- rep(lwd, length.out = length(cv))

 # border and colors 
  Col <- createcolors(facets, border, Col)
      
 # update polygons.
  poly <- 
      list(x       = cbind(poly$x, PolyX),
           y       = cbind(poly$y, PolyY),               
           z       = cbind(poly$z, PolyZ),               
           col     = c(poly$col, Col$facet),
           border  = c(poly$border, Col$border),
           lty     = c(poly$lty, Lty),
           lwd     = c(poly$lwd, Lwd),
           isimg   = c(poly$isimg, rep(0, length.out = length(cv))),
           img     = poly$img
        )

  return(poly)
} 

## =============================================================================
## Adding polygons from images or persps
## =============================================================================

addimg <- function(poly, x, y, z, colvar = z, plist,
                   col = NULL, NAcol = "white", breaks = NULL,
                   border = NA, facets = TRUE, lwd = 1, lty = 1,
                   resfac = 1, clim = NULL,   
                   ltheta = -135, lphi = 0, shade = NA, 
                   lighting = FALSE, alpha = NULL, ...)  {

  if (missing(poly) | is.null(poly) | length(poly) == 0) 
    poly <- list(x = NULL, y = NULL, col = NULL, border = NULL, 
                 lwd = NULL, lty = NULL, proj = NULL, isimg = NULL,
                 img = list())                    

  else if (class(poly) != "poly")
    stop ("'poly' not of correct type for addimg")

  dot <- splitdotpersp(list(), 
                       bty = NULL, lighting, x, y, z, plist = plist, 
                       ltheta = ltheta, lphi = lphi, shade = shade, breaks = breaks)

  if (! is.null(plist$xs)) {
    dot$shade$xs <- plist$scalefac$x
    dot$shade$ys <- plist$scalefac$y
    dot$shade$zs <- plist$scalefac$z
  }  

  if (! is.matrix(z))
    stop("'z' should be a matrix")

  if (! is.null(colvar)) {
    if (! is.matrix(colvar))
      stop("'colvar' should be a matrix or NULL")

    if (any (dim(z) != dim(colvar)) )
      stop("dimension of 'colvar' not equal to dimension of 'z'")

    if (is.null(col))
      col <- jet.col(100)
  }

  if (is.null(col))
    col <- "grey"
  
 # if x and y are a vector: check resfac and convert to matrix
  X <- x
  Y <- y

  if (is.vector(X))  {
    if (! is.vector(Y))
      stop("'y' should be a vector if 'x' is one")

    if (any(resfac != 1)) {   # change resolution
      res <- changeres(resfac, X, Y, z, colvar)
      X <- res$x
      Y <- res$y
      z <- res$z
      colvar <- res$colvar
    }

    XY <- mesh(X, Y)
    x <- XY$x
    y <- XY$y
  } else { 
    x <- X
    y <- Y
  }
     
 # check class and dimensionality
  if (! is.matrix(x))
    stop("'x' should be a matrix or vector")
  if (! is.matrix(y))
    stop("'y' should be a matrix or vector")

  DD <- dim(x)
  if (any (DD != dim(y)) )
    stop("dimension of 'x' not equal to dimension of 'y'")
  if (any (DD != dim(z)) )
    stop("dimension of 'y' not equal to dimension of 'z'")

 # create polygons - unsorted                  
  ix <- rep(1:nrow(x), ncol(x))
  iy <- as.vector(matrix(nrow = nrow(x), ncol = ncol(x), 
                     data = 1:ncol(x), byrow =TRUE))

  Poly <- createpoly(x, y, z, ix, iy) 

 # depth view of the points 
  Proj   <- project(colMeans(Poly$X, na.rm = TRUE), 
                    colMeans(Poly$Y, na.rm = TRUE), 
                    colMeans(Poly$Z, na.rm = TRUE), plist)

 # colvar is converted to colors.
  if (! is.null(colvar)) {
    if (is.null(clim))
      clim <- range(colvar, na.rm = TRUE)     

 # Colors for values = NA 
    if (any (is.na(colvar)) & ! is.null(NAcol) ) {
      CC <- checkcolors(colvar, col, NAcol, clim)
      clim   <- CC$lim
      col    <- CC$col
      colvar <- CC$colvar
    }
    if (is.null(breaks)) {
      crange <- diff(clim)
      N      <- length(col) -1
      Col <- col[1 + trunc((colvar - clim[1])/crange*N)]
    } else {
      zi <- .bincode(colvar, breaks, TRUE, TRUE)
      Col <- col[zi]
    }
  } else 
    Col <- rep(col, length.out = length(x))
  Col[is.na(Col)] <- NAcol

  if (! dot$shade$type == "none") 
    Col <- facetcols (x, y, z, Col, dot$shade)

  if (! is.null(alpha)) 
    Col <- setalpha(Col, alpha)

  imgCol <- createcolors(facets, border, Col)
  
  if (is.null(alpha)) alpha <- NA
  Alpha <- alpha
  alpha <- rep(alpha, length.out = length(x))
  
 # update polygons.
  numimg <- length(poly$img)    
  poly <- 
    list(x      = cbind(poly$x, Poly$X),
         y      = cbind(poly$y, Poly$Y),               
         z      = cbind(poly$z, Poly$Z),               
         col    = c(poly$col, imgCol$facet),
         border = c(poly$border, imgCol$border),
         lwd    = c(poly$lwd, rep(lwd , length.out = length(x))),
         lty    = c(poly$lty, rep(lty , length.out = length(x))),
         proj   = c(poly$proj, Proj),
         isimg  = c(poly$isimg, rep(1, length.out = length(x))),
         alpha  = c(poly$alpha, alpha),
         img    = poly$img)
  if (numimg == 0)
    poly$img <- list()
  sl <- list (ix = ix, iy = iy, Proj = Proj, list = 1:length(ix))
  poly$img[[numimg+1]] <- list(x = X, y = Y, z = z, 
               col = matrix(nrow = nrow(colvar), 
               ncol = ncol(colvar), data = Col), NAcol = NAcol,
               breaks = breaks, sl = sl,
               facets = facets, border = border, lwd = lwd, lty = lty,
               alpha = Alpha, mapped = TRUE)      

  class(poly) <- "poly"
  
  return(poly)
}


## =============================================================================
## Adds lines to segments
## =============================================================================

addlines <- function(segm, x, y, z, plist,
                     col, NAcol = "white", 
                     lwd = 1, lty = 1, 
                     ltheta = -135, lphi = 0, shade = NA, 
                     lighting = FALSE, ignorez = TRUE, ...)
                      {

  if (missing(segm) | is.null(segm)) 
    segm <- list(x.from = NULL, x.to = NULL, 
                 y.from = NULL, y.to = NULL, 
                 z.from = NULL, z.to = NULL,
                 col = NULL, lwd = NULL, lty = NULL, proj = NULL)                    
  else if (class(segm) != "segments")
    stop ("'segm' not of correct type, 'segments', for addlines")
                             
  dot <- splitdotpersp(list(), 
                       bty = NULL, lighting, x, y, z, plist = plist, 
                       ltheta = ltheta, lphi = lphi, shade = shade, breaks = NULL)
  if (is.null(col))
    col <- "black"
  
 # check class and dimensionality
  if (! is.vector(x))
    stop("'x' should be a vector")
  if (! is.vector(y))
    stop("'y' should be a vector")
  if (! is.vector(z))
    stop("'z' should be a vector")

  len <- length(x)
  if (len != length(y) )
    stop("dimension of 'x' not equal to dimension of 'y'")
  if (len != length(z) )
    stop("dimension of 'x' not equal to dimension of 'z'")

 # depth view of the midpoints of each segment 

  Proj   <- project(0.5*(x[-1] +x[-len]), 0.5*(y[-1]+y[-len]), 
    0.5*(z[-1]+z[-len]), plist, ignorez = ignorez)
  Col <- rep(col, length.out = len-1)

  if (! dot$shade$type == "none") 
    Col <- facetcols (x, y, z, Col, dot$shade)

  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = len-1)
  
 # update and return segments.
  segm <- list(
       x.from = c(segm$x.from, x[-len]),
       x.to   = c(segm$x.to,   x[-1]),
       y.from = c(segm$y.from, y[-len]),                                  
       y.to   = c(segm$y.to,   y[-1]),                                  
       z.from = c(segm$z.from, z[-len]),                                  
       z.to   = c(segm$z.to,   z[-1]),                                  
       col    = c(segm$col, Col),
       lwd    = c(segm$lwd, rep(lwd , length.out = len-1)),
       lty    = c(segm$lty, rep(lty , length.out = len-1)),
       alpha  = c(segm$alpha, alpha),
       proj   = c(segm$proj, Proj))
  class(segm) <- "segments"
  return(segm)
}

## =============================================================================

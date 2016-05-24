## =============================================================================
## 3-D histograms, x, y matrix or vector; z = matrix
## =============================================================================

hist3D <- function(x = seq(0, 1, length.out = nrow(z)),
                   y = seq(0, 1, length.out = ncol(z)),
                   z, ..., 
                   colvar = z, phi = 40, theta = 40,
                   col = NULL, NAcol = "white", breaks = NULL,
                   border = NA, facets = TRUE, colkey = NULL,
                   image = FALSE, contour = FALSE, panel.first = NULL,
                   clim = NULL, clab = NULL, bty = "b",
                   lighting = FALSE, shade = NA, ltheta = -135, lphi = 0,
                   space = 0, opaque.top = FALSE,
                   add = FALSE, plot = TRUE) {

  if (! is.matrix(z))
    stop("'z'  should be a matrix") 

  if (! is.vector(x) & length(dim(x)) != 1)
    stop("'x' should be a vector")
    
  if (! is.vector(y) & length(dim(y)) != 1)
    stop("'y' should be a vector")

  if (length(x) != nrow(z))
    stop("'x' should be of length = nrow(z)")
  
  if (length(y) != ncol(z))
    stop("'y' should be of length = ncol(z)")

  if (any(space > 0.9))
    stop(" 'space' too large, should be smaller or equal to 0.9")

  else if (any(space < 0.0))
    stop(" 'space' cannot be smaller than 0")

  space <- rep(space, length.out = 2) / 2

  plist <- initplist(add)
  ll <- list(...)
  dot <- splitdotpersp(ll, bty, lighting, 
    extendvec(x), extendvec(y), z, plist = plist, shade, lphi, ltheta,
    breaks = breaks)
  
 # swap if decreasing
  if (length(x) > 0 & all(diff(x) < 0)) {     
    x <- rev(x)
    if (is.null(dot$persp$xlim)) 
      dot$persp$xlim <- range(x)
      else if (diff(dot$persp$ylim) < 0)
        stop("'persp' expects increasing ylim")  
      else if (diff(dot$persp$xlim) < 0)
        stop("'persp' expects increasing xlim")  

    if (ispresent(colvar)) 
      colvar <- colvar[nrow(colvar):1, ]
    z <- z[nrow(z):1, ]
  }
 
  if (length(y) > 1 & all(diff(y) < 0)) {     
    y <- rev(y)
    if (is.null(dot$persp$ylim)) 
      dot$persp$ylim <- range(y)
    else if (diff(dot$persp$ylim) < 0)
      stop("'persp' expects increasing ylim")  
    if (ispresent(colvar)) 
      colvar <- colvar[, (ncol(colvar):1)]
    z <- z[, (ncol(z):1)]
  }

  image   <- check.args(image)
  contour <- check.args(contour)
  if (contour$add) 
    cv <- colvar

  if (is.null(col) & is.null(breaks))
   col <- jet.col(100)
  else if (is.null(col))
   col <- jet.col(length(breaks)-1)
  breaks <- check.breaks(breaks,col)

  CC <- check.colvar.2(colvar, z, col, clim, dot$alpha) 
  topcol <- (! is.null(dot$alpha) & opaque.top)
  if (topcol)  
    col2 <- check.colvar.2(colvar, z, col, clim, NULL)$col
  
  colvar <- CC$colvar
  col <- CC$col

  if (ispresent(colvar)) {
    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
       
    if (dot$clog) {                    
      colvar <- log(colvar)
      clim   <- log(clim)
    }
  
    iscolkey <- is.colkey(colkey, col) 
    if (iscolkey)
      colkey <- check.colkey(colkey)

 # colors
    crange <- diff(clim)
    N      <- length(col) - 1
    cmin   <- clim[1]

    if (is.null(breaks))
     Cols <- matrix(nrow = nrow(z),
       col[1 + trunc((colvar - cmin)/crange*N)])
    else
     Cols <- matrix(nrow = nrow(z),
       col[.bincode(colvar, breaks, TRUE, TRUE)])
    Cols [is.na(Cols)] <- NAcol
    if (topcol) {
      if (is.null(breaks))
       Col2 <- matrix(nrow = nrow(z),
         col2[1 + trunc((colvar - cmin)/crange*N)])
      else
       Col2 <- matrix(nrow = nrow(z),
         col2[.bincode(colvar, breaks, TRUE, TRUE)])
      Col2 [is.na(Col2)] <- NAcol
     }  
  } else { 
    iscolkey <- FALSE
    Cols <- rep(col , length.out = length(z))
    if (topcol) 
      Col2 <- rep(col2, length.out = length(z))
  }
  
 # mapping from centre to interfaces
  extend <- function(x) {
   # This does exceed the x- y boundaries
    N <- length(x)          
    x <- c(x[1] - (x[2]-x[1]), x, x[N] + (x[N]-x[(N-1)]))  
    ii <- 2:length(x)
    0.5*(x[ii] + x[ii-1])           
  }

 # expand x and y to grid 
  XYmesh <- mesh(x, y)

  if (length (x) > 1)
    xx <- extend(x)
  else 
    xx <- dot$persp$xlim
      
  if (length (y) > 1)
    yy <- extend(y)
  else 
    yy <- dot$persp$ylim
 
  if (is.null(plist)) {
    do.call("perspbox", c(alist(range(xx), range(yy), dot$persp$zlim, 
                   phi = phi, theta = theta, plot = plot, 
                   colkey = colkey, col = col), dot$persp))
    plist <- getplist()
  }
    
  if (! is.null(panel.first))
    panel.first(plist$mat)

 # viewing order
  ind <- expand.sort(1:length(XYmesh$x), dim(XYmesh$x)) 
  ix <- ind$x; iy <- ind$y
  
 # The colors
  Col <- createcolors(facets, border, Cols)
  if (topcol) 
    Col2 <- createcolors(facets, border, Col2)

 # the polygons:
 #     2
 #  3     4
 #     1
 # 5 = top polygon
 
 # shading?
  isshade <- dot$shade$type != "none"
  islight <- FALSE

  if (isshade) {

    to5facets <- function(val)
      as.vector(matrix(ncol = 5, 
        data = rep(val, length(ix)), byrow = TRUE))   

    light   <- setuplight(dot$shade$lphi, dot$shade$ltheta) 

  # the normals are known for rectangles:
    norms <- matrix(nrow = 5, byrow = TRUE, data =
       c(0, -1, 0,   0, 1, 0,   -1, 0, 0,   1, 0, 0,   0, 0, 1)  )
                                                                                                                   
    if (dot$shade$type == "shade") {
      shade <- abs(dot$shade$shade) 
      Sum <- 0.5*(norms[,1]*light[1] + norms[,2]*light[2] + norms[,3]*light[3] +1)
      Shade <- to5facets (Sum^shade)
    } else  {
      isshade <- FALSE
      islight <- TRUE
      Normals <- list(u = to5facets(norms[, 1]), 
                      v = to5facets(norms[, 2]), 
                      w = to5facets(norms[, 3]))
    }
  }  
  PolyX <- NULL
  PolyY <- NULL
  PolyZ <- NULL
  COL <- BORD <- NULL
  
  dx <- diff(xx)*space[1]
  dy <- diff(yy)*space[2]

 # basal and top points of the column; x and y positions
  z.k    <- rep(min(dot$persp$zlim), length(z))
  z.kp1  <- as.vector(z)
  x.i    <- xx[ix  ]+dx[ix]
  x.ip1  <- xx[ix+1]-dx[ix]
  y.j    <- yy[iy  ]+dy[iy]
  y.jp1  <- yy[iy+1]-dy[iy]

 # the facet coordinates
  PolyX <- cbind(rbind(x.i  , x.ip1, x.ip1, x.i  ),
                 rbind(x.ip1, x.i  , x.i  , x.ip1),
                 rbind(x.i  , x.i  , x.i  , x.i  ),
                 rbind(x.ip1, x.ip1, x.ip1, x.ip1),
                 rbind(x.i  , x.ip1, x.ip1, x.i  ))
  PolyY <- cbind(rbind(y.j  , y.j  , y.j  , y.j  ),
                 rbind(y.jp1, y.jp1, y.jp1, y.jp1),
                 rbind(y.j  , y.jp1, y.jp1, y.j  ),
                 rbind(y.j  , y.jp1, y.jp1, y.j  ),
                 rbind(y.j  , y.j  , y.jp1, y.jp1))

  PolyZ <- cbind(rbind(z.k  , z.k  , z.kp1, z.kp1),
                 rbind(z.k  , z.k  , z.kp1, z.kp1),
                 rbind(z.k  , z.k  , z.kp1, z.kp1),
                 rbind(z.k  , z.k  , z.kp1, z.kp1),
                 rbind(z.kp1, z.kp1, z.kp1, z.kp1))

 # facet colors  
  if (!topcol) 
    COL   <- rep(Col$facet , 5)#[i,j] 
  else
    COL   <- c(rep(Col$facet , 4), Col2$facet)#[i,j] 
    
  if (!topcol) 
    BORD  <- rep(Col$border, 5)#[i,j] 
  else
    BORD  <- c(rep(Col$border, 4), Col2$border)#[i,j] 
  
  if (isshade) {
     if (facets) {
       RGB  <- t(col2rgb(COL)) * Shade / 255
       COL  <- rgb(RGB)
       
     if (! is.null(dot$alpha)) {
         if (!topcol) 
           COL <- setalpha(COL, dot$alpha)
         else {
           ialph <- 1:(4*length(Col$facet))
           COL[ialph] <- setalpha(COL[ialph], dot$alpha)
         }    
     }           
             
     }
     if (! is.na(border)){
       RGB  <- t(col2rgb(BORD)) * Shade / 255
       BORD <- rgb(RGB) 
       if (! is.null(dot$alpha)) BORD <- setalpha(BORD, dot$alpha)
     } 
  } else if (islight) {
     if (facets) 
       COL  <- facetcols.shadelight (light, Normals, COL,  dot$shade)
     if (! is.na(border))
       BORD <- facetcols.shadelight (light, Normals, BORD, dot$shade)
  }

  PolyX <- rbind(PolyX, NA)   
  PolyY <- rbind(PolyY, NA)   
  PolyZ <- rbind(PolyZ, NA)   

  Proj <- project(colMeans(PolyX, na.rm = TRUE),
                  colMeans(PolyY, na.rm = TRUE), 
                  colMeans(PolyZ, na.rm = TRUE), plist, TRUE)

  lwd <- dot$points$lwd
  if (is.null(lwd)) 
    lwd <- 1

  lty <- dot$points$lty
  if (is.null(lty)) 
    lty <- 1
  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = ncol(PolyX))
    
  Poly <- list(x      = PolyX,
               y      = PolyY,                                  
               z      = PolyZ,
               col    = COL,
               border = BORD,
               lwd    = rep(lwd , length.out = ncol(PolyX)),
               lty    = rep(lty , length.out = ncol(PolyX)),
               isimg  = rep(0, length.out = ncol(PolyX)),
               alpha  = alpha, 
               proj   = Proj)
  class(Poly) <- "poly"

  if (image$add) 
    Poly <- XYimage (Poly, image, x, y, z, plist, col, breaks = breaks)

  if (contour$add) 
    segm <- contourfunc(contour, x, y, z, plist, cv, clim)
  else
    segm <- NULL

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
      dot$clog, type = "hist3D", breaks = breaks)

  plist <- plot.struct.3D(plist, poly = Poly, segm = segm, plot = plot)  

  setplist(plist)
  invisible(plist$mat)

}


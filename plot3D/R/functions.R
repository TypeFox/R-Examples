## =============================================================================
## =============================================================================
## General utility functions
## =============================================================================
## =============================================================================

## =============================================================================
## Extend/average to embrace points, but not exceeding boundaries
## =============================================================================

extendvec <- function(x) {
  ll <- length(x)  
  c(x[1] + 0.5*(x[1] - x[2]), x, x[ll] + 0.5*(x[ll] - x[ll-1]))           
}


extend <- function(x, na.rm = TRUE) {
  if (na.rm & any(is.na(x)))
    return(extend.na(x))
  x <- rbind(x[1,], x, x[nrow(x),])  
  x <- cbind(x[,1], x, x[,ncol(x)])
  ii <- 2:nrow(x)
  jj <- 2:ncol(x)           
  0.25*(x[ii, jj] + x[ii-1, jj] + x[ii, jj-1] + x[ii-1, jj-1])           
}

# extending, but ignoring the NAs
extend.na <- function(x) {
  x <- rbind(x[1,], x, x[nrow(x),])  
  x <- cbind(x[,1], x, x[,ncol(x)])
  nisna <- !is.na(x)
  x[!nisna] <- 0
  ii <- 2:nrow(x)
  jj <- 2:ncol(x)           
  Sum <- (x[ii, jj] + x[ii-1, jj] + x[ii, jj-1] + x[ii-1, jj-1])
  Count <- nisna[ii, jj] + nisna[ii-1, jj] + nisna[ii, jj-1] + nisna[ii-1, jj-1]
  Res <- Sum/Count
  Res[is.nan(Res)] <- NA
  return(Res)
}

meangrid.na <- function(x) {
  nisna <- !is.na(x)
  x[!nisna] <- 0
  nr <- nrow(x)
  nc <- ncol(x)
  Sum <-  x[-1, -1] + x[-1, -nc] + x[-nr, -1] + x[-nr, -nc]
  Count <- nisna[-1, -1] + nisna[-1, -nc] + nisna[-nr, -1] + nisna[-nr, -nc]
  Res <- Sum/Count
  Res[is.nan(Res)] <- NA
  return(Res)
  
}
meangrid <- function(x, na.rm = FALSE) {
  if (na.rm & any(is.na(x)))
    return(meangrid.na(x))
  else
    return(0.25*( x[-1,       -1] + x[-1,       -ncol(x)] 
                + x[-nrow(x), -1] + x[-nrow(x), -ncol(x)]))
}


## =============================================================================
## something can be toggled off by putting it = FALSE, NA, NULL
## =============================================================================

ispresent <- function(var) {
  if (length(var) > 1) 
    return(TRUE)

  ispresent <- TRUE
  if (is.null(var))
    ispresent <- FALSE
  else if (is.na(var))
    ispresent <- FALSE
  else if (length(var) == 1)
    if (is.logical(var)) 
      if (!var)
        ispresent <- FALSE
  return(ispresent)
}

## =============================================================================
## Segment function that allows for multiple lwd  - not necessary? -not used
## =============================================================================

Segments <- function (x0, y0, x1, y1, col, lty, lwd, ...) {
  nD <- length(x0)

  doseg <- function(i, ...) {
    lines (x = c(x0[i], x1[i]), y = c(y0[i], y1[i]), 
           col = col[i], lty = lty[i], lwd = lwd[i], ...)
  }
  lapply (1:nD, FUN = doseg)

}

## =============================================================================
## List to draw legend/contour 
## =============================================================================

check.args <- function(ll) {
  
  addit <- ! is.null(ll)
 
  if (length(ll) == 0) 
    addit <- FALSE
  
  else if (is.logical(ll[[1]]))
    if (length(ll[[1]]) == 1)      
      if (ll[[1]] == FALSE) 
        addit <- FALSE
  
  side <- "zmin"
  
  if (addit) {     # should have at least side argument
    if (is.list(ll)) {
      if (!is.null(ll$side))
        side <- ll$side
      ll$side <- NULL  
    } else ll <- list() 
  }
  
  list(add = addit, side = side, args = ll)
}

## =============================================================================
## Segments to plot contours on bottom or top panel
## =============================================================================
    
contourfunc <- function(contour, x, y, z, plist, cv = NULL, 
  clim = range(cv), dDepth = NULL, addbox = TRUE) { 

  if (is.null(dDepth))
    dDepth <- contour$args$dDepth
  if (is.null(dDepth))
    dDepth <- 1e-1
  levels <- contour$args$levels
  nlevels <- contour$args$nlevels
  if (is.null (nlevels  ))
    nlevels <- 10
  if (! is.null(contour$args$x)) 
    x <- contour$args$x
  if (! is.null(contour$args$y)) 
    y <- contour$args$y
  if (! is.null(contour$args$z)) 
    z <- contour$args$z
  contour$args$x <- contour$args$y <- contour$args$z <- NULL
  if (is.null(cv))
    cv <- z
  contour$args$lighting <- contour$args$shade <- NULL
  contour$args$levels <- contour$args$nlevels <- NULL
  
  if (! is.null(levels))
    line.list <- 
      contourLines(x, y, cv, levels = levels)
  else
    line.list <-contourLines(x, y, cv, nlevels = nlevels)

  col <- contour$args$col

  if (is.null(col))
    col <- "black"
    
  if (length(col) == 1)
    col <- c(col, col)  

  crange <- diff(clim)
  N      <- length(col) -1

  getcol <- function(v) 
     col[1 + trunc((v - clim[1])/crange*N)]

  contour$args$col <- contour$args$dDepth <- NULL
  segm <- NULL
  
  for (side in contour$side) {  
  
    if (side == "zmin")
      zz <- min(plist$zlim)
    else if (side == "zmax")
      zz <- max(plist$zlim)
    else if (side == "z")
      zz <- NULL
    else if (!is.numeric(as.numeric(side)))
      stop ("cannot add contour on side ", side)
    else
      zz <- as.numeric(side)

    if (side == "z")  {     # contour *on* the persp plot

      Nx <- length(x)
      Ny <- length(y)
      dx <- c(diff(x), 1)  # 1= for last value
      dy <- c(diff(y), 1)

      for (i in 1:length(line.list)) {     
        xto <- line.list[[i]]$x
        yto <- line.list[[i]]$y

   # find embracing values : first interval
        ix <- findInterval(xto, x)
        iy <- findInterval(yto, y)

   # next interval
        ixp1 <- pmin(ix+1, Nx)
        iyp1 <- pmin(iy+1, Ny)

   # interpolation factor
        xfac <- (xto-x[ix])/dx[ix]
        yfac <- (yto-y[iy])/dy[iy]

   # interpolate
        zz <- (1-yfac)*((1-xfac)*z[cbind(ix,iy)  ]+xfac*z[cbind(ixp1,iy)]) +
                  yfac*((1-xfac)*z[cbind(ix,iyp1)]+xfac*z[cbind(ixp1,iyp1)])
        Col <- getcol(line.list[[i]]$level)

        segm <- do.call("addlines", c(alist(segm, line.list[[i]]$x, line.list[[i]]$y, 
          z = as.vector(zz), col = Col, plist = plist), contour$args))
      }      

    } else
      for (i in 1:length(line.list))      
        segm <- do.call("addlines", c(alist(segm, line.list[[i]]$x, 
          line.list[[i]]$y, z = rep(zz, length(line.list[[i]]$x)), 
          col = getcol(line.list[[i]]$level), plist = plist), contour$args))

     if (side != "z" & addbox)      
      segm <- addlines(segm, x = c(x[1], x[length(x)],x[length(x)], x[1], x[1]),
               y = c(y[1], y[1], y[length(y)],y[length(y)], y[1]), 
               z = rep(zz, length.out = 5), col = "black", plist = plist)
  }
  segm$proj <- segm$proj + dDepth*plist$persp$expand  # put it on foreground...

  return(segm)
}
          
## =============================================================================
## polygons if image is to be drawn on bottom or top panel
## =============================================================================

XYimage <- function(poly, image, x, y, z,  plist, col, breaks) {

  if (is.null(image$args$col))
    image$args$col <- col

  if (! is.null(image$args$x)) 
    x <- image$args$x
  if (! is.null(image$args$y)) 
    y <- image$args$y
  if (! is.null(image$args$z)) 
    z <- image$args$z
  image$args$x <- image$args$y <- image$args$z <- NULL
  
  lwd <- image$args$lwd
  if (is.null(lwd)) 
    lwd <- 1

  lty <- image$args$lty
  if (is.null(lty)) 
    lty <- 1

  image$args$lwd <- image$args$lty <- NULL

  for (side in image$side) {  

    if (side == "zmin")
      zz <- plist$zlim[1]

    else if (side == "zmax")
      zz <- plist$zlim[2]

    else if (!is.numeric(as.numeric(side)))
      stop ("cannot add image on side ", side)

    else
      zz <- as.numeric(side)

    zmat <- matrix(nrow = length(x), ncol = length(y), data = zz)           
    poly <- do.call("addimg", c(alist(poly, x, y, z = zmat, 
        colvar = z, breaks = breaks, plist = plist, lwd = lwd, lty = lty), image$args))
        
  }                                                
  return(poly)
}

## =============================================================================
## Functions that account for occurrence of decreasing values...
## =============================================================================

FindInterval <- function(x, vec, ...) {

  if (all(diff(vec) < 0)) { 
    vec <- rev(vec)
    res <- c(length(vec):1) [findInterval(x, vec, ...)]-1
  } else 
    res <- findInterval(x, vec, ...)
  res [ res == 0] <- 1
  res
}

## .. and of NAs
Approx <- function(x, y, ...) {

  if (all(is.na(x))) 
    return(list(y = rep(NA, length = length(y)), x = y))

  if (diff(range(x, na.rm = TRUE)) == 0)
    warning("Warning in approx: all 'x' values are the same")

  if (any(is.na(c(x, y)))) {
    ii <- unique(c(which(is.na(x)), which(is.na(y))))
    x <- x[-ii]
    y <- y[-ii]
  } 
  approx(x, y, ...)
}

## =============================================================================
## =============================================================================
## expands a sorted index list to matrix rows and columns 
## =============================================================================
## =============================================================================

expand.sort <- function(sortlist, Dim) {

 # can be improved!
  II <- matrix(nrow = Dim[1], ncol = Dim[2], data = 1:Dim[1])
  ix <- II[sortlist]    # indices to x-values in sorted list
  II <- matrix(nrow = Dim[1], ncol = Dim[2], data = 1:Dim[2], byrow =TRUE)
  iy <- II[sortlist]    # indices to y-values in sorted list
  list(x = ix, y = iy)
}

## =============================================================================
## =============================================================================
## Projection depth 
## =============================================================================
## =============================================================================

project <- function (x, y, z, plist, ignorez = TRUE) {
 
  if (is.null(plist))
    stop("a 3D plot has not yet been created")
  if (ignorez) z <- 0           # it is done like this in persp (plot.c)
  TV     <- cbind(as.vector(x), as.vector(y), as.vector(z), 1) %*% plist$mat
  return (-TV[ ,4])
}
    
## =============================================================================
## =============================================================================
## wrapper over trans3d that returns arrays with original dimensions
## =============================================================================
## =============================================================================

trans3D <- function(x, y, z, pmat) {

 # trans3d expects vectors as input
  XX <- trans3d(x = as.vector(x), y = as.vector(y), z = as.vector(z), 
                pmat = pmat)
   
 # convert x- and y to an array with original dimensions
  if (! is.vector(x)) {
    XX$x <- array(dim = dim(x), data = XX$x)              
    XX$y <- array(dim = dim(x), data = XX$y)              
  } 
  return(XX)
}
    
## =============================================================================
## =============================================================================
## Function to change the resolution of matrices - 
## =============================================================================
## =============================================================================

changeres <- function(resfac, x, y, z, colvar = NULL, na.rm = FALSE) { 

  if (is.matrix(x)) 
    return(changeres_mat(resfac, x, y, z, colvar, na.rm))
  else 
    return(changeres_vec(resfac, x, y, z, colvar, na.rm))
}

changeres_vec <- function(resfac, x, y, z, colvar = NULL, na.rm = FALSE) { 
    
  resfac <- abs(rep(resfac, length.out = 2))
  diffx <- diff(x)
  diffy <- diff(y)
  XX <- x
  YY <- y 
  RX <- 1/resfac[1]
  RY <- 1/resfac[2]
  if (resfac[1] > 1)
    for (i in 1: (resfac[1]-1))  
      XX <- c(XX, x[-length(x)] + diffx * i*RX)
  else if (resfac[1] < 0.99) 
    XX <- x[as.integer(seq(1, nrow(z), length.out = nrow(z)*resfac[1]))]
  if (resfac[2] > 1)
    for (i in 1: (resfac[2]-1))  
      YY <- c(YY, y[-length(y)] + diffy * i*RY)
  else if (resfac[2] < 0.99) 
    YY <- y[as.integer(seq(1, ncol(z), length.out = ncol(z)*resfac[2]))]

  XX <- unique(sort(XX))
  YY <- unique(sort(YY))
  z <- remapxy(z, x = x, y = y, xto = XX, yto = YY, na.rm)
  
  if (! is.null(colvar))  
    colvar <- remapxy(colvar, x = x, y = y, xto = XX, yto = YY, na.rm)
  
  list(x = XX, y = YY, z = z, colvar = colvar)
}

changeres_mat <- function(resfac, x, y, z, colvar = NULL, na.rm = FALSE) { 

  resfac <- abs(rep(resfac, length.out = 2))
  xx <- 1:nrow(z)
  yy <- 1:ncol(z)
  XX <- changeres(resfac, xx, yy, x)$z
  YY <- changeres(resfac, xx, yy, y)$z
  RR  <- changeres(resfac, xx, yy, z, colvar)
  
  list(x = XX, y = YY, z = RR$z, colvar = RR$colvar)
}

## =============================================================================
## Maps a matrix 'z' from (x, y) to (xto, yto) by linear 2-D interpolation
## =============================================================================

remapxy <- function(z, x, y, xto, yto, na.rm = FALSE) {     

  if (na.rm & any(is.na(z)))
    return(remapxyNA (z, x = x, y = y, xto = xto, yto = yto))

# a simple function with linear interpolation - only for x and y a vector
  Nx <- length(x)
  Ny <- length(y)
    
  dx  <- c(diff(x), 1)  # 1= for last value
  dy  <- c(diff(y), 1)

 # find embracing values : first interval
  ix <- FindInterval(xto, x)
  iy <- FindInterval(yto, y)

 # interpolation factor
  xfac <- (xto-x[ix])/dx[ix]
  yfac <- (yto-y[iy])/dy[iy]

 # expand for all combinations..
  gg <- expand.grid(ix, iy)
  ix <- gg[,1]
  iy <- gg[,2]
    
 # next interval
  ixp1 <- pmin(ix+1, Nx)
  iyp1 <- pmin(iy+1, Ny)

  gg <- expand.grid(xfac, yfac)
  xfac <- gg[,1]
  yfac <- gg[,2]

 # interpolate
  M <-(1-yfac)*((1-xfac)*z[cbind(ix,iy)  ]+xfac*z[cbind(ixp1,iy)]) +
          yfac*((1-xfac)*z[cbind(ix,iyp1)]+xfac*z[cbind(ixp1,iyp1)])
  
  return (matrix(nrow = length(xto), ncol = length(yto), data = M))
}

remapxyNA <- function(z, x, y, xto, yto) {
    Nx <- length(x)
    Ny <- length(y)
    dx <- c(diff(x), 1)
    dy <- c(diff(y), 1)
    ix <- FindInterval(xto, x)
    iy <- FindInterval(yto, y)
    xfac <- (xto - x[ix])/dx[ix]
    yfac <- (yto - y[iy])/dy[iy]
    gg <- expand.grid(ix, iy)
    ix <- gg[, 1]
    iy <- gg[, 2]
    ixp1 <- pmin(ix + 1, Nx)
    iyp1 <- pmin(iy + 1, Ny)
    gg <- expand.grid(xfac, yfac)
    xfac <- gg[, 1]
    yfac <- gg[, 2]

    f <- zz <- matrix(nrow = length(xfac), ncol = 4)
    zz[ , 1] <- z[cbind(ix, iy)]
    f[ , 1] <- (1 - yfac) *(1 - xfac)

    zz[,2]  <- z[cbind(ix, iyp1)]
    f[,2] <- yfac * (1 - xfac)

    f[,3] <- (1 - yfac) * xfac
    zz[,3]  <- z[cbind(ixp1, iy)]
 
    f[,4] <- yfac * xfac
    zz[,4] <- z[cbind(ixp1, iyp1)]

    naii <- is.na(zz)
    f[naii] <- 0
    zz[naii] <- 0
     
    rows <- rowSums(f)
    f <- f/rowSums(f)
    
    M <- f[,1] * zz[,1] + f[,2] * zz[,2] + f[,3] * zz[,3] + f[,4] * zz[,4]
    return(matrix(nrow = length(xto), ncol = length(yto), data = M))
}

## =============================================================================
## =============================================================================
## Colors
## =============================================================================
## =============================================================================

## =============================================================================
## Generates color vector based on variable values
## =============================================================================

variablecol <- function(colvar, col, NAcol, clim, breaks) {
 
  if (is.null(breaks)) {
    ncol <- length(col)
  
    colvar[colvar < min(clim)] <- NA
    colvar[colvar > max(clim)] <- NA
    rn <- clim[2] - clim[1]
    ifelse (rn != 0, Col <- col[1 + trunc((colvar - clim[1])/rn *
      (ncol - 1)+1e-15)], Col <- rep(col[1], ncol))               # + tiny: since R 3.2.2
  } else {
      zi <- .bincode(colvar, breaks, TRUE, TRUE)
      Col <- col[zi]
  }
  Col[is.na(Col)] <- NAcol
  return(Col)
}

## =============================================================================
## Check/extend colors for NAs and adapt range
## =============================================================================

checkcolors <- function(colvar, col, NAcol, lim) {

  colvar[colvar < min(lim)] <- NA
  colvar[colvar > max(lim)] <- NA

  if (length(col) == 1) 
    col <- c(col, col)
  N <- length(col) - 1
  col <- c(NAcol, col)
  rr <- diff(lim)
  
  if (rr == 0) 
    rr <- lim[1] *1e-3
  colvar[is.na(colvar)] <- lim[1] - 1/N*rr
  lim [1] <- lim[1] - 1/N*rr
  
  list(col = col, colvar = colvar, lim = lim)
}

## =============================================================================
## Colors for facets and border
## =============================================================================

createcolors <- function(isfacets, border, Cols) {

  isnaborder <- FALSE
  if (length(border) > 0)
    isnaborder <- is.na(border)

  if (is.na(isfacets)) {
    if (isnaborder)        
      border <- Cols
    Cols[] <- NA

  } else if (isfacets) {    # facets added
    if (isnaborder) {
      border <- Cols
      border[] <- NA 
    }
  } else {
    if (is.na(border))      # no facets
      border <- Cols
    Cols[] <- "white"
  }
  if (length(border ) == 1) {
    bb <- border
    border <- Cols
    border[] <- bb
  }
  
  list(border = border, facet = Cols)
}
  
## =============================================================================
## check dimensionality of colvar and colors - for inttype = 2
## =============================================================================

check.colvar.2 <- function(colvar, z, col, clim, alpha) {
  
  iscolvar <- ispresent(colvar)

  if (iscolvar) {
    
    if (any (dim(colvar) - dim(z)) != 0)
      stop("dimension of 'colvar' should be equal to dimension of 'z'")    
    
    if (! is.null(clim)) {
      if (length(clim) != 2)
        stop("'clim' should be a two-valued vector with the ranges of 'colvar'")
      colvar[colvar < min(clim)] <- NA
      colvar[colvar > max(clim)] <- NA
    }  
            
   # check colors  
    if (length(col) == 1) 
      if (is.na(col)) 
        col <- NULL

    if (is.null(col))
      col <- jet.col(100)

    if (length(col) == 1) {
      col <- c(col, col)
    }  

  } else { 
     if (is.null(col))
      col <- rep("grey", 2)
    else
      col <- rep(col[1], 2)  # take first color
  }
  if (! is.null(alpha))
     col <- setalpha(col, alpha)
  
  list(colvar = colvar, col = col)     
}

## =============================================================================
## check dimensionality of colvar and colors
## =============================================================================

check.colvar.persp <- function(colvar, z, col, inttype, clim, alpha) {
  
  iscolvar <- ispresent(colvar)

  if (iscolvar) {
    
    if (inttype == 2 & any (dim(colvar) - dim(z)) != 0)
      stop("dimension of 'colvar' should be equal to dimension of 'z'")    
    
    else if (inttype != 2){
      if (all (dim(colvar) - dim(z)) == 0)
        colvar <- meangrid(colvar, inttype == 3)  # averages of colvar
      else if (any (dim(colvar) - dim(z) != -1))
        stop("dimension of 'colvar' should be equal to dimension of 'z' or have one row and one column less")    
    }
    
    if (! is.null(clim)) {
      if (length(clim) != 2)
        stop("'clim' should be a two-valued vector with the ranges of 'colvar'")
      colvar[colvar < min(clim)] <- NA
      colvar[colvar > max(clim)] <- NA
    }  
            
   # check colors  
    if (length(col) == 1) 
      if (is.na(col)) 
        col <- NULL

    if (is.null(col))
      col <- jet.col(100)

    if (length(col) == 1) {
      col <- c(col, col)
    }  

  } else { 
     if (is.null(col))
      col <- rep("grey", 2)
    else
      col <- rep(col[1], 2)  # take first color
  }
  if (! is.null(alpha))
     col <- setalpha(col, alpha)
  
  list(colvar = colvar, col = col)     
}

## =============================================================================
## ranges, scaling factors
## =============================================================================

setlim <- function(xlim, ylim, zlim, scale, expand) {
  
  if (is.null(scale))
    scale <- TRUE
  if (is.null(expand))
    expand <- 1

  xs <- 0.5 *abs(diff(xlim))
  ys <- 0.5 *abs(diff(ylim))
  zs <- 0.5 *abs(diff(zlim))
                      
  xc <- 0.5 *sum(xlim)
  yc <- 0.5 *sum(ylim)
  zc <- 0.5 *sum(zlim)

  if (! scale) {
    ss <- max (xs, ys, zs)
    xs <- ss; ys <- ss; zs <- ss
  }
  xs <- ifelse (xs == 0, 1, 1 / xs)
  ys <- ifelse (ys == 0, 1, 1 / ys)
  zs <- ifelse (zs == 0, 1, expand / zs)
  list(x = xs, y = ys, z = zs, xc = xc, yc = yc, zc = zc, expand = expand)
}  

## =============================================================================
## Shade and lighting-related parameters
## =============================================================================

check.shade <- function(shadedots, lighting) {

  if (is.null(shadedots$lphi))
    shadedots$lphi <- 0
  if (is.na(shadedots$lphi))
    shadedots$lphi <- 0
  if (is.null(shadedots$ltheta))
    shadedots$ltheta <- -135
  if (is.na(shadedots$ltheta))
    shadedots$ltheta <- -135
   
  shadedots$type <- "none"
  if (! is.null(lighting)) {         
    if (is.character(lighting))
      shadedots$type <- "light"
    else if (is.logical(lighting)) {
      if (lighting)
        shadedots$type <- "light"
    } else if (is.list(lighting)) {
      if (!is.null(lighting$type)) 
        shadedots$type <- lighting$type
      else   
        shadedots$type <- "light"
      lighting$type <- NULL
      shadedots <- c(shadedots, lighting)
    }  
  }  
   
  if (! is.null(shadedots$shade))
    if (! is.na(shadedots$shade) & shadedots$type == "none") # lighting overrules shade
      shadedots$type <- "shade"     
  if (is.null(shadedots$shade)) 
    shadedots$shade <- NA
  return (shadedots)
}  

## =============================================================================
## Split dots in part related to persp/other + set clog, labels
## =============================================================================

splitdotpersp <- function(dots, bty = "b", lighting = NULL, 
   x = NULL, y = NULL, z = NULL, plist = NULL,
   shade = NA, lphi = 0, ltheta = -135, breaks) {

  dots$bty <- bty

  namespersp <- c("xlim", "ylim", "zlim", "xlab", "ylab", "zlab",
        "main", "sub", "r", "d", "scale", "expand","box", "axes", 
        "nticks", "ticktype", "col.ticks", "lwd.ticks", "bty", 
        "cex.axis", "col.axis", "font.axis", "xaxs", "yaxs", 
        "col.panel", "lwd.panel", "col.grid", "lwd.grid",
        "cex.lab", "col.lab", "font.lab",   
        "cex.main", "col.main", "font.main", "alpha")
            
  setlim <- c(!is.null(dots$xlim), !is.null(dots$ylim), ! is.null(dots$zlim))

  # log of color variable
  clog <- dots$clog
  if (is.null(clog)) { 
    clog <- FALSE
    if (! is.null(dots[["log"]])) {
      if (length(grep("c", dots[["log"]])) > 0) {
        dots[["log"]] <- gsub("c", "", dots[["log"]])
        if (dots[["log"]] == "")
          dots[["log"]] <- NULL
        clog <- TRUE
      } 
    } 
  }
  if (clog & ! is.null(breaks)) {
    warning("cannot combine log = 'c' and 'breaks' - removing log = 'c'")
    clog <- FALSE
  }

 # labels        
  if (is.null(dots$xlab)) 
    dots$xlab <- "x"
  if (is.null(dots$ylab))     
    dots$ylab <- "y"
  if (is.null(dots$zlab)) 
    dots$zlab <- "z"

 # ranges  
  if (! is.null(plist)) {
      dots$xlim <- plist$xlim
      dots$ylim <- plist$ylim
      dots$zlim <- plist$zlim                    
      scalefac <- plist$scalefac
  } else if (! is.null(x)) {
    if (is.null(dots$xlim)) 
      dots$xlim <- range(x, na.rm = TRUE)
    if (is.null(dots$ylim)) 
      dots$ylim <- range(y, na.rm = TRUE)
    if (is.null(dots$zlim)) 
      dots$zlim <- range(z, na.rm = TRUE)
    scalefac <- setlim (dots$xlim, dots$ylim, dots$zlim, dots$scale, dots$expand)
  } else 
    scalefac <- list(x = NULL, y = NULL, z = NULL, xc = NULL, yc = NULL, zc = NULL)
  
 # shade and lighting parameters
  shadedots <- list(ltheta = ltheta, lphi = lphi, shade = shade)
  shadedots <- check.shade(shadedots, lighting)
    
  Persp <- c(dots[ names(dots) %in% namespersp], clog = clog, setlim = setlim)  
  if (! is.null(dots$alpha)) {
    if (! is.numeric(dots$alpha))
      stop("'alpha' should be numeric") 
    if (length(dots$alpha) > 1)
      stop("'alpha' should be one number") 
    if (dots$alpha < 0 | dots$alpha > 1)
      stop("'alpha' should be a number inbetween 0 and 1") 
  }  
  list(persp = Persp,
       points = dots[!names(dots) %in% c(namespersp, "clog", "alpha")],
       shade = c(shadedots, xs = scalefac$x, ys = scalefac$y, zs = scalefac$z,
         alpha = dots$alpha), 
       clog = clog, alpha = dots$alpha)
}      

## =============================================================================
## Expanding arguments in dots  (...)
## =============================================================================

repdots <- function(dots, n) {
  if (is.function(dots)) { 
    return(dots)
  } else 
    return(rep(dots, length.out = n))
}

expandsortdots <- function(dot, sortlist) {
  ls <- length(sortlist)
  expsort <- function(dot) {
    if (length(dot) > 1)
      rep(dot, length.out = ls) [sortlist]
    else
      dot
  }      
  lapply(dot, expsort) 
}

setdots <- function(dots, n) 
  lapply(dots, repdots, n)

## =============================================================================
## Extracting element 'index' from dots  (...)
## =============================================================================

extractdots <- function(dots, index) {
  
  ret <- lapply(dots, "[", index)
  ret <- lapply(ret, unlist) # flatten list
  return(ret)
}
## =============================================================================
## set mfrow and ask
## =============================================================================

setplotpar <- function(ldots, nv, ask) {

  nmdots <- names(ldots) 

  # nv = number of variables to plot
  if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
    nc <- min(ceiling(sqrt(nv)), 3)
    nr <- min(ceiling(nv/nc), 3)
    mfrow <- c(nr, nc)
  } else if ("mfcol" %in% nmdots)
    mfrow <- rev(ldots$mfcol)
  else 
    mfrow <- ldots$mfrow

  if (! is.null(mfrow))  
    mf <- par(mfrow = mfrow)

  ## interactively wait if there are remaining figures
  if (is.null(ask))
    ask <- prod(par("mfrow")) < nv && dev.interactive()

  return(ask)
}

## =============================================================================
## Split plotting parameters in general (main) and point parameters
## =============================================================================

splitpardots <- function(dots) {

  clog <- dots$clog
  if (is.null(clog)) { 
    clog <- FALSE
    if (! is.null(dots$log)) {
      if (length(grep("c", dots[["log"]])) > 0) {
        dots[["log"]] <- gsub("c", "", dots[["log"]])
        if (dots[["log"]] == "")
          dots[["log"]] <- NULL
        clog <- TRUE
      } 
    } 
  }
  
  nmdots <- names(dots)

  # plotting parameters : split in plot parameters and point parameters
  plotnames <- c("xlab", "ylab", "zlab", "xlim", "ylim", "zlim", 
                 "main", "sub", "log", "asp", "bty", 
                 "xaxs", "yaxs", "xaxt", "yaxt", "breaks",
                 "ann", "axes", "frame.plot", "panel.first", "panel.last",
                 "cex.lab", "col.lab", "font.lab", "las", "tck", "tcl", "mgp", 
                 "cex.axis", "col.axis", "font.axis", 
                 "cex.main", "col.main", "font.main")

  # plot.default parameters
  ii <- names(dots) %in% plotnames
  dotmain <- dots[ii]

  # point parameters
  ip <- !names(dots) %in% c(plotnames, "add", "clog", "alpha")
  dotpoints <- dots[ip]
  # alpha
  if (! is.null(dots$alpha)) {
    if (! is.numeric(dots$alpha))
      stop("'alpha' should be numeric") 
    if (length(dots$alpha) > 1)
      stop("'alpha' should be one number") 
    if (dots$alpha < 0 | dots$alpha > 1)
      stop("'alpha' should be a number inbetween 0 and 1") 
  }  
  
  list (points = dotpoints, main = dotmain, add = dots$add, 
    clog = clog, alpha = dots$alpha)

}


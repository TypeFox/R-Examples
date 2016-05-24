## =============================================================================
## Mapping, sections, selection
## =============================================================================
remap <- function(var, ...) UseMethod ("remap")
extract <- function(var, ...) UseMethod ("extract")
changeres <- function(var, ...) UseMethod("changeres")

## =============================================================================
## Increase resolution from a matrix by a factor - x can be a matrix or a vector
## =============================================================================
changeres.matrix <- function(var, x, y, resfac, na.rm = TRUE, ...) { 

  resfac <- abs(rep(resfac, length.out = 2))
  if (is.matrix(x) | is.matrix(y)) 
    return(changeres_xmat(resfac, x, y, var, na.rm))
  else 
    return(changeres_xvec(resfac, x, y, var, na.rm))
}

changeres.array <- function(var, x, y, z, resfac, na.rm = TRUE, ...) { 

  if (! is.vector(x) | ! is.vector(y) | ! is.vector(z))
    stop ("'x', 'y', and 'z' should be a vector")

  resfac <- abs(rep(resfac, length.out = 3))
  xto <- changeresvec(x, resfac[1])
  yto <- changeresvec(y, resfac[2])
  zto <- changeresvec(z, resfac[3])
  
  return(remap.array (var, x, y, z, xto, yto, zto, na.rm, ...))
}

changeres_xvec <- function(resfac, x, y, z, na.rm) { 
    
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
  if (na.rm & any(is.na(z)))
    z <- remapxyNA(z, x, y, XX, YY)
  else
    z <- remapxy(z, x, y, XX, YY)
  
  list(var = z, x = XX, y = YY)
}

changeres_xmat <- function(resfac, x, y, z, na.rm) { 

  xx <- 1:nrow(z)
  yy <- 1:ncol(z)
  XX <- changeres_xvec(resfac, xx, yy, x, na.rm)$var 
  YY <- changeres_xvec(resfac, xx, yy, y, na.rm)$var
  ZZ <- changeres_xvec(resfac, xx, yy, z, na.rm)$var
  
  list(var = ZZ, x = XX, y = YY)
}


## =============================================================================
## change resolution to arbitrary set of x, y (z) values
## =============================================================================

changeresvec <- function(x, resfac) {
  resfac[is.na(resfac)] <- 1

  diffx <- diff(x) 
  xto <- x
  RX <- 1/resfac
  Nx <- length(x)
  if (resfac > 1)
    for (i in 1: (resfac-1))  
      xto <- c(xto, x[-length(x)] + diffx * i*RX)

  else if (resfac < 0.99) 
    xto <- x[as.integer(seq(1, Nx, length.out = Nx*resfac))]

  xto <- sort(xto)
  return(xto)
}

# same as from package plot3D....
remapxy <- function (z, x, y, xto, yto) {
    Nx <- length(x)    
    Ny <- length(y)
    
    dx <- c(diff(x), 1)     # 1 for last value
    dy <- c(diff(y), 1)

 # find embracing values : first interval
    ix <- FindInterval(xto, x)
    iy <- FindInterval(yto, y)
 # interpolation facotr
    xfac <- (xto - x[ix])/dx[ix]
    yfac <- (yto - y[iy])/dy[iy]
 # expand for all combinations..
    gg <- expand.grid(ix, iy)
    ix <- gg[, 1]
    iy <- gg[, 2]
 # next inetrval
    ixp1 <- pmin(ix + 1, Nx)
    iyp1 <- pmin(iy + 1, Ny)
    gg <- expand.grid(xfac, yfac)
    xfac <- gg[, 1]
    yfac <- gg[, 2]
 # interpolate
    M <- (1 - yfac) * ((1 - xfac) * z[cbind(ix, iy)] + xfac * 
        z[cbind(ixp1, iy)]) + yfac * ((1 - xfac) * z[cbind(ix, 
        iyp1)] + xfac * z[cbind(ixp1, iyp1)])
    return(matrix(nrow = length(xto), ncol = length(yto), data = M))
}

## =============================================================================

remapxyNA <- function(z, x = x, y = y, xto = x, yto = y) {
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
    zz[, 1] <- z[cbind(ix, iy)]
    f[, 1]  <- (1 - yfac) *(1 - xfac)

    zz[, 2] <- z[cbind(ix, iyp1)]
    f[, 2]  <- yfac * (1 - xfac)

    f[, 3]  <- (1 - yfac) * xfac
    zz[, 3] <- z[cbind(ixp1, iy)]
 
    f[, 4]  <- yfac * xfac
    zz[, 4] <- z[cbind(ixp1, iyp1)]

    naii <- is.na(zz)
    f[naii] <- 0
    zz[naii] <- 0
     
    rows <- rowSums(f)
    f <- f/rowSums(f)
    
    M <- rowSums(f * zz)
    return(matrix(nrow = length(xto), ncol = length(yto), data = M))
}

## =============================================================================
## 2-D mapping, x, y a vector, var a matrix
## =============================================================================

remap.matrix <- function(var, x, y, xto = NULL, yto = NULL, na.rm = TRUE,...) {

  if (is.array(x)) {
    if (length(dim(x)) !=  1)
      stop ("'x' should be a vector or array of dimension 1")
    x <- as.vector(x)  
  }
  if (is.array(y)) {
    if (length(dim(y)) !=  1)
      stop ("'y' should be a vector or array of dimension 1")
    y <- as.vector(y)  
  }
  if (! is.vector(x) | ! is.vector(y))
    stop ("'x' and 'y' should be a vector")
    
  if (is.array(var)) {
    if (length(dim(var)) !=  2)
      stop ("'var' should be a matrix or array of dimension 2")
    var <- as.matrix(var)  
  }

  Nx <- length(x)
  Ny <- length(y)

  DD <- dim(var)
  
  if (is.null(xto))
    xto <- seq(min(x), max(x), length.out = Nx)
  if (is.null(yto))
    yto <- seq(min(y), max(y), length.out = Ny)

  if (min(xto) < min(x) | max(xto) > max(x))
    stop("'x' should embrace 'xto'")
  if (min(yto) < min(y) | max(yto) > max(y))
    stop("'y' should embrace 'yto'")

  Du <- dim(var)

  if (length(Du) != 2)
    stop ("'var' should be a matrix")
  if (Du[1] != Nx)
    stop("'var' and 'x' not compatible: 1st dimension not equal to length (x)")
  if (Du[2] != Ny)
    stop("'var' and 'y' not compatible: 2nd dimension not equal to length (y)")

  if (na.rm & any(is.na(var)))
    M <- remapxyNA(var, x, y, xto, yto)
  else
    M <- remapxy(var, x, y, xto, yto)

  list (var = M, x = as.vector(xto), y = as.vector(yto))
}

## =============================================================================
## THREE DIMENSIONS
## =============================================================================

remapxyz <- function(var, x, y, z, xto, yto, zto) {
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  dx  <- c(diff(x), 1)  # 1= for last value
  dy  <- c(diff(y), 1)
  dz  <- c(diff(z), 1)

 # find embracing values : first interval
  ix <- FindInterval(xto, x)
  iy <- FindInterval(yto, y)
  iz <- FindInterval(zto, z)

 # interpolation factors
  xfac <- (xto-x[ix])/dx[ix]
  yfac <- (yto-y[iy])/dy[iy]
  zfac <- (zto-z[iz])/dz[iz]

 # expand for all combinations..
  gg <- expand.grid(ix,iy,iz)
  ix <- gg[,1]
  iy <- gg[,2]
  iz <- gg[,3]
    
  # next interval
  ixp1 <- pmin(ix+1, Nx)
  iyp1 <- pmin(iy+1, Ny)
  izp1 <- pmin(iz+1, Nz)

  gg <- expand.grid(xfac,yfac,zfac)
  xfac <- gg[,1]
  yfac <- gg[,2]
  zfac <- gg[,3]

  # interpolate
   MM <-   (1 - zfac) *
      ((1-yfac)*((1-xfac)*var[cbind(ix,iy,iz)  ]+xfac*var[cbind(ixp1,iy,iz)]) +
          yfac *((1-xfac)*var[cbind(ix,iyp1,iz)]+xfac*var[cbind(ixp1,iyp1,iz)])
      ) + zfac* (
      (1-yfac)*((1-xfac)*var[cbind(ix,iy,izp1)  ]+xfac*var[cbind(ixp1,iy,izp1)]) +
          yfac*((1-xfac)*var[cbind(ix,iyp1,izp1)]+xfac*var[cbind(ixp1,iyp1,izp1)]))
   return(array(dim = c(length(xto), length(yto), length(zto)), data = MM))
}

## =============================================================================

remapxyzNA <- function(var, x, y, z, xto, yto, zto) {
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)

  dx  <- c(diff(x), 1)  # 1= for last value
  dy  <- c(diff(y), 1)
  dz  <- c(diff(z), 1)

 # find embracing values : first interval
  ix <- FindInterval(xto, x)
  iy <- FindInterval(yto, y)
  iz <- FindInterval(zto, z)

 # interpolation factors
  xfac <- (xto-x[ix])/dx[ix]
  yfac <- (yto-y[iy])/dy[iy]
  zfac <- (zto-z[iz])/dz[iz]

 # expand for all combinations..
  gg <- expand.grid(ix,iy,iz)
  ix <- gg[,1]
  iy <- gg[,2]
  iz <- gg[,3]
    
  # next interval
  ixp1 <- pmin(ix+1, Nx)
  iyp1 <- pmin(iy+1, Ny)
  izp1 <- pmin(iz+1, Nz)

  gg <- expand.grid(xfac,yfac,zfac)
  xfac <- gg[,1]
  yfac <- gg[,2]
  zfac <- gg[,3]

  # interpolate
  f <- zz <- matrix(nrow = length(xfac), ncol = 8)
  zz[ , 1] <- var[cbind(ix, iy, iz)]
  f[ , 1]  <- (1-zfac) * (1-yfac) * (1-xfac)
  zz[ , 2] <- var[cbind(ixp1, iy, iz)]
  f[ , 2]  <- (1-zfac) * (1-yfac) * (xfac) 
  zz[ , 3] <- var[cbind(ix, iyp1, iz)]
  f[ , 3]  <- (1-zfac) * (yfac) * (1-xfac)
  zz[ , 4] <- var[cbind(ixp1, iyp1, iz)]
  f[ , 4]  <- (1-zfac) * (yfac) * (xfac) 
  zz[ , 5] <- var[cbind(ix, iy, izp1)]
  f[ , 5]  <- (zfac) * (1-yfac) * (1-xfac)
  zz[ , 6] <- var[cbind(ixp1, iy, izp1)]
  f[ , 6]  <- (zfac) * (1-yfac) * (xfac) 
  zz[ , 7] <- var[cbind(ix, iyp1, izp1)]
  f[ , 7]  <- (zfac) * (yfac) * (1-xfac)
  zz[ , 8] <- var[cbind(ixp1, iyp1, izp1)]
  f[ , 8]  <- (zfac) * (yfac) * (xfac) 

  naii <- is.na(zz)
  f[naii] <- 0
  zz[naii] <- 0
     
  rows <- rowSums(f)
  f <- f/rowSums(f)
  
  MM <- colSums(f * zz)

   return(array(dim = c(length(xto), length(yto), length(zto)), data = MM))
}

## =============================================================================
## 3-D mapping, x, y, z a vector; var a 3-D array.
## =============================================================================

remap.array <- function(var, x, y, z, xto = NULL, yto = NULL, 
  zto = NULL, na.rm = TRUE, ...) {

  if (is.array(x)) {
    if (length(dim(x)) !=  1)
      stop ("'x' should be a vector or array of dimension 1")
    x <- as.vector(x)  
  }
  if (is.array(y)) {
    if (length(dim(y)) !=  1)
      stop ("'y' should be a vector or array of dimension 1")
    y <- as.vector(y)  
  }
  if (is.array(z)) {
    if (length(dim(z)) !=  1)
      stop ("'z' should be a vector or array of dimension 1")
    z <- as.vector(z)  
  }

  if (! is.vector(x) | ! is.vector(y) | ! is.vector(z))
    stop ("'x', 'y', and 'z' should be a vector")

  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  DD <- dim(var)
  
  if (is.null(xto))
    xto <- seq(min(x), max(x), length.out = Nx)
  if (is.null(yto))
    yto <- seq(min(y), max(y), length.out = Ny)
  if (is.null(zto))
    zto <- seq(min(z), max(z), length.out = Nz)

  if (min(xto) < min(x) | max(xto) > max(x))
    stop("'x' should embrace 'xto'")
  if (min(yto) < min(y) | max(yto) > max(y))
    stop("'y' should embrace 'yto'")
  if (min(zto) < min(z) | max(zto) > max(z))
    stop("'z' should embrace 'zto'")

  Du <- dim(var)

  if (length(Du) != 3)
    stop ("'var' should be an array of dim 3")
  if (Du[1] != Nx)
    stop("'var' and 'x' not compatible: 1st dimension not equal to length (x)")
  if (Du[2] != Ny)
    stop("'var' and 'y' not compatible: 2nd dimension not equal to length (y)")
  if (Du[3] != Nz)
    stop("'var' and 'z' not compatible: 3rd dimension not equal to length (z)")

  if (na.rm & any(is.na(z)))
    M <- remapxyzNA(var, x, y, z, xto, yto, zto)
  else
    M <- remapxyz(var, x, y, z, xto, yto, zto)
  
  list (var = M, x = as.vector(xto), y = as.vector(yto), z = as.vector(zto))
}

## =============================================================================
## =============================================================================
## Takes a selection across a matrix 'var' from (x, y) to 
## cbind(xto, yto) by linear 2-D interpolation
## x, y: vector
## =============================================================================
## =============================================================================

extract.matrix <- function(var, x, y, xyto, ...) {

  Nx <- length(x)
  Ny <- length(y)

  if (ncol(xyto) != 2)
    stop("'xyto' should be a two-columned matrix with x and y values")
  
  xto <- xyto[ ,1]  
  yto <- xyto[ ,2]
    
  if (min(xto) < min(x) | max(xto) > max(x)) 
    stop("'x' should embrace elements in first column of 'xyto'")
  if (min(yto) < min(y) | max(yto) > max(y)) 
    stop("'y' should embrace elements in second column of 'xyto'")

  dx  <- c(diff(x), 1)  # 1= for last value
  dy  <- c(diff(y), 1)

  Du <- dim(var)
  if (length(Du) != 2)
    stop ("'var' should be a matrix")
  if (Du[1] != Nx)
    stop("'var' and 'x' not compatible: 1st dimension not equal to length (x)")
  if (Du[2] != Ny)
    stop("'var' and 'y' not compatible: 2nd dimension not equal to length (y)")

 # find embracing values : first interval
  ix <- FindInterval(xto, x )
  iy <- FindInterval(yto, y)

 # next interval
  ixp1 <- pmin(ix+1, Nx)
  iyp1 <- pmin(iy+1, Ny)

 # interpolation factor
  xfac <- (xto-x[ix])/dx[ix]
  yfac <- (yto-y[iy])/dy[iy]

 # interpolate
  MM <- (1-yfac)*((1-xfac)*var[cbind(ix,iy)  ]+xfac*var[cbind(ixp1,iy)]) +
            yfac*((1-xfac)*var[cbind(ix,iyp1)]+xfac*var[cbind(ixp1,iyp1)])

  colnames(xyto) <- c("x", "y")
  list (var = MM, xy = xyto)             
}

## =============================================================================
## =============================================================================
## Takes a selection across an array 'var' from (x, y, z) to 
## cbind(xto, yto, zto) by linear interpolation
## x, y and z: vector
## =============================================================================
## =============================================================================

extract.array <- function(var, x, y, z, xyzto, ...) {

  if (! is.vector(x) | ! is.vector(y) | ! is.vector(z))
    stop ("'x', 'y', and 'z' should be a vector")

  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  
  if (ncol(xyzto) != 3)
    stop("'xyzto' should be a three-columned matrix with x, y, and z- values")
  
  xto <- xyzto[ ,1]  
  yto <- xyzto[ ,2]
  zto <- xyzto[ ,3]
    
  if (min(xto) < min(x) | max(xto) > max(x)) 
    stop("'x' should embrace elements in first column of 'xyzto'")
  if (min(yto) < min(y) | max(yto) > max(y)) 
    stop("'y' should embrace elements in second column of 'xyzto'")
  if (min(zto) < min(z) | max(zto) > max(z)) 
    stop("'z' should embrace elements in third column of 'xyzto'")

  dx  <- c(diff(x), 1)  # 1= for last value
  dy  <- c(diff(y), 1)
  dz  <- c(diff(z), 1)

  Du <- dim(var)
  if (length(Du) != 3)
    stop ("'var' should be an array of dim 3")
  if (Du[1] != Nx)
    stop("'var' and 'x' not compatible: 1st dimension not equal to length (x)")
  if (Du[2] != Ny)
    stop("'var' and 'y' not compatible: 2nd dimension not equal to length (y)")
  if (Du[3] != Nz)
    stop("'var' and 'z' not compatible: 3rd dimension not equal to length (z)")

 # find embracing values : first interval
  ix <- FindInterval(xto, x )
  iy <- FindInterval(yto, y)
  iz <- FindInterval(zto, z)
  
 # next interval
  ixp1 <- pmin(ix+1, Nx)
  iyp1 <- pmin(iy+1, Ny)
  izp1 <- pmin(iz+1, Nz)
  
 # interpolation factor
  xfac <- (xto-x[ix])/dx[ix]
  yfac <- (yto-y[iy])/dy[iy]
  zfac <- (zto-z[iz])/dz[iz]
  
 # interpolate
   MM <-   (1 - zfac) *
      ((1-yfac)*((1-xfac)*var[cbind(ix,iy,iz)  ]+xfac*var[cbind(ixp1,iy,iz)]) +
          yfac *((1-xfac)*var[cbind(ix,iyp1,iz)]+xfac*var[cbind(ixp1,iyp1,iz)])
      ) + zfac* (
      (1-yfac)*((1-xfac)*var[cbind(ix,iy,izp1)  ]+xfac*var[cbind(ixp1,iy,izp1)]) +
          yfac*((1-xfac)*var[cbind(ix,iyp1,izp1)]+xfac*var[cbind(ixp1,iyp1,izp1)]))

  colnames(xyzto) <- c("x", "y", "z")
  list (var = MM, xyz = xyzto)             
}

## =============================================================================
## =============================================================================
## Takes a transect across an array 'var' 
## =============================================================================
## =============================================================================

transect <- function(var, x, y, z, to, margin = "xy", ...) {

  if (! is.vector(x) | ! is.vector(y) | ! is.vector(z))
    stop ("'x', 'y', and 'z' should be a vector")

  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)

  Du <- dim(var)  
  if (length(Du) != 3)
    stop ("'var' should be an array")
  if (Du[1] != Nx)
    stop("'var' and 'x' not compatible: 1st dimension not equal to length (x)")
  if (Du[2] != Ny)
    stop("'var' and 'y' not compatible: 2nd dimension not equal to length (y)")
  if (Du[3] != Nz)
    stop("'var' and 'z' not compatible: 3rd dimension not equal to length (z)")

  if (ncol(to) != 2)
    stop("'to' should be a two-columned matrix with ", margin, " values")
  
#  if (!margin %in% c("xy", "yx", "xz", "zx", "yz", "zy"))
#    stop ("elements in 'margin' should be either 'x', 'y', 'z', and there should be two of them")
  if (!margin %in% c("xy", "xz", "yz"))
    stop ("'margin' should be one of 'xy', 'xz', 'yz'")

# test for validity of values
  testmargin <- function(i, margin, to) {
    if (margin == "x") 
      if (min(to) < min(x) | max(to) > max(x)) 
        stop("'x' should embrace elements in column ", i, " of 'to'")
    else if (margin == "y") 
      if (min(to) < min(y) | max(to) > max(y)) 
        stop("'y' should embrace elements in column ", i, " of 'to'")
    else if (margin == "z") 
      if (min(to) < min(z) | max(to) > max(z)) 
        stop("'z' should embrace elements in column ", i, " of 'to'")
  }
  testmargin(1, substr(margin, 1, 1), to[ ,1])
  testmargin(2, substr(margin, 2, 2), to[ ,2])
  
  along <- c("z", "y", "x")[match(margin, c("xy", "xz", "yz"))]
  
  if (along == "x") {
    xyzto <- cbind(expand.grid(x, to[,1]), to[,2]) 
    result <- extract(var, x, y, z, xyzto)$var
    result <- list (var = t(matrix (ncol = nrow(to), data = result)),
                    rows = to, cols = x)   
    colnames(result$rows) <- c("y","z")

  }  else if (along == "y") {
#    stop ("not yet implemented for margin xz")
    Var <- aperm(var, c(1, 3, 2))
    xyzto <- cbind(to[,1], expand.grid(to[,2], y))
    result <- extract(Var, x, z, y, xyzto)$var
    result <- list (var = matrix (nrow = nrow(to), data = result),
                    rows = to, cols = y)   
    colnames(result$rows) <- c("x","z")

  } else if (along == "z") {
    xyzto <- cbind(to[,1], expand.grid(to[,2], z))
    result <- extract(var, x, y, z, xyzto)$var
    result <- list (var = matrix (nrow = nrow(to), data = result),
                    rows = to, cols = z)   
    colnames(result$rows) <- c("x","y")
  }
  result
}
## =============================================================================
## =============================================================================
## Maps or extracts from variables defined in sigma coordinates 
## =============================================================================
## =============================================================================

mapsigma <- function(var, ...) UseMethod ("mapsigma")

## =============================================================================

mapsigma.matrix <- function (var = NULL, 
                             sigma,      
                             signr = 2,
                             x = NULL,
                             depth = NULL,
                             numdepth = NULL, 
                             xto = NULL,
                             resfac = 1, ...) {
  if (is.null(signr)) 
    signr <- 2

  if (is.null(var)) 
    var <- matrix(nrow = nrow(sigma), ncol = ncol(sigma), data = 1)

  if (any (dim(var) - dim(sigma) != 0))
    stop ("'sigma' and 'var' not of same dimension")
    
  DD <- dim(var)
  Dxy <- DD[-signr]
  if (is.null(x))
    x <- seq(0, 1, length.out = Dxy[1])

  if (length(x) != Dxy[1])
    stop("dimensions of 'var' and 'x' not compatible: ", Dxy[1], " not = ", length(x))

  yto <- y <- seq(0, 1, length.out = DD[signr] )
  if (length(y) != DD[signr])
    stop("dimensions of 'var' and 'y' not compatible: ", DD[signr], " not = ", length(y))

  resfac <- abs(rep(resfac, length.out = 2))
  changeres <- FALSE
  if (any(resfac > 1)) {
    if (is.null(xto))
      xto <- changeresvec(1:Dxy, resfac[1])
    yto <- changeresvec(1:DD[signr], resfac[2])
    changeres <- TRUE
  } else 
    if (! is.null(xto) )
      changeres <- TRUE

  if (changeres) {
  dots <- list(...)
  dots$x <- dots$y <- NULL
    var <- do.call("remap", c(alist(var, x = 1:Dxy, y = 1:DD[signr], xto = xto, yto = yto), dots))$var
    sigma <- do.call("remap", c(alist(sigma, x = 1:Dxy, y = 1:DD[signr], xto = xto, yto = yto), ...))$var
    xy <- (1:2)[-signr]
    x <- changeresvec(x, resfac[xy])
  }
  
  Nr <- nrow(var)
  Nc <- ncol(var)
    
  if (is.null(depth))
    depth <- seq(min(sigma, na.rm = TRUE), max(sigma, na.rm = TRUE), 
               length.out = ifelse(is.null(numdepth),ncol(sigma),numdepth))

  if (nrow(sigma) != Nr | ncol(sigma) != Nc)
    stop ("'sigma' should be of same dimension as matrix 'var'")

  if (signr == 2) { 
    Mnew <- matrix(nrow = Nr, ncol = length(depth), data =NA)
  
    for (i in 1:Nr) 
      if (any(! is.na(sigma[i, ])))
        Mnew[i,] <- Approx(x = sigma[i,], y = var[i,], xout = depth)$y 
  } else {
    Mnew <- matrix(nrow = length(depth), ncol = Nc, data =NA)
  
    for (i in 1:Nc) 
      if (any(! is.na(sigma[ ,i])))
        Mnew[ ,i] <- Approx(x = sigma[ ,i], y = var[,i], xout = depth)$y 
  }   
  list (var = Mnew, depth = depth, x = x)
}

## =============================================================================

mapsigma.array <- function (var = NULL, 
                            sigma,
                            signr = 3,
                            x = NULL,
                            y = NULL, 
                            depth = NULL,
                            numdepth = NULL, 
                            xto = NULL,
                            yto = NULL,
                            resfac = 1, ...) {
  if (is.null(signr)) 
    signr <- 3

  if (is.null(var)) 
    var <- array(dim = dim(sigma), data = 1)

  DD <- dim(var)
  Ds <- dim(sigma)
    
  if (any (DD - Ds != 0))
    stop ("'sigma' should be of same dimension as 'var'")
    
  lenD <- length(DD)
  if (lenD > 3)
    stop ("dimension of 'sigma' or 'var' can not be > 3")
  
  Dxy <- DD[-signr]
  if (is.null(x))
    x <- seq(0, 1, length.out = Dxy[1])

  if (is.null(y))
    y <- seq(0, 1, length.out = Dxy[2])
  
  zto <- z <- seq(0, 1, length.out = DD[signr])

  resfac <- abs(rep(resfac, length.out = 3))
  changeres <- FALSE
  if (any(resfac > 1)) {
    if (is.null(xto))
      xto <- changeresvec(1:Dxy[1], resfac[1])
    if (is.null(yto))
      yto <- changeresvec(1:Dxy[2], resfac[2])
    zto <- changeresvec(1:DD[signr], resfac[3])
    changeres <- TRUE
  } else 
    if (! is.null(xto) | ! is.null(yto) )
      changeres <- TRUE

  if (changeres) { 
    var <- remap(var, x = 1:Dxy[1], y = 1:Dxy[2], z = 1:DD[signr], 
               xto = xto, yto = yto, zto = zto, ...)$var
    sigma <- remap(sigma, x = 1:Dxy[1], y = 1:Dxy[2], z = 1:DD[signr], 
               xto = xto, yto = yto, zto = zto, ...)$var
    xy <- (1:3)[-signr]
    x <- changeresvec(x, resfac[xy[1]])
    y <- changeresvec(y, resfac[xy[2]])
    
  }
  DD <- dim(var)
  Dxy <- DD[-signr]

  if (is.null(depth))
    depth <- seq(min(sigma, na.rm = TRUE), max(sigma, na.rm = TRUE), 
                length.out = ifelse(is.null(numdepth), dim(sigma)[signr], numdepth))

  
  if (signr == 3) {
    Mnew <- array(dim = c( Dxy[1:2], length(depth)), data = NA)
    for (i in 1:Dxy[1])
      for (j in 1:Dxy[2]) {
        if (any(! is.na(sigma[i,j,])))
          Mnew[i,j,] <- Approx(x = sigma[i,j,], y = var[i,j,], xout = depth)$y 
      }
  } else if (signr == 2) {
    Mnew <- array(dim = c(Dxy[1], length(depth), Dxy[2]), data = NA)
    for (i in 1:Dxy[1])
      for (j in 1:Dxy[2]) {
        if (any(! is.na(sigma[i,,j])))
          Mnew[i,,j] <- Approx(x = sigma[i,,j], y = var[i,,j], xout = depth)$y 
      }
  } else  if (signr == 1) {
    Mnew <- array(dim = c( length(depth), Dxy[1:2]), data = NA)
    for (i in 1:Dxy[1])
      for (j in 1:Dxy[2]) {
        if (any(! is.na(sigma[,i,j])))
          Mnew[,i,j] <- Approx(x = sigma[,i,j], y = var[,i,j], xout = depth)$y 
      }
  }

  list (var = Mnew, depth = depth, x = x, y = y)
} 
## =============================================================================
## Transect in sigma coordinates
## =============================================================================

transectsigma  <- function (var = NULL, 
                            sigma,
                            x = NULL,
                            y = NULL,
                            to, 
                            depth = NULL,
                            numdepth = NULL,
                            resfac = 1, ...) {

  if (is.null(var)) 
    var <- array(dim = dim(sigma), data = 1)

  DD <- dim(var)
  Ds <- dim(sigma)
    
  if (any (DD - Ds != 0))
    stop ("'sigma' should be of same dimension as 'var'")
    
  lenD <- length(DD)
  if (lenD != 3)
    stop ("dimension of 'sigma' or 'var' should be 3")
    
  z <-  seq(0, 1, length.out = DD[3])


  tranVar <- transect(var = var, x = as.vector(x),
     y = as.vector(y), z = z, to = to)

  tranSigma <- transect(var = sigma, x = as.vector(x),
    y = as.vector(y), z = z, to = to)

  VarSigma <- mapsigma(tranVar$var, x = to[,1], y = to[,2],
     sigma = tranSigma$var, depth = depth, numdepth = numdepth, resfac = resfac, ...)
  
  return(VarSigma)
  
}
## =============================================================================
## Creates a rectangular grid from 2 or 3 vectors (x, y, z)
## =============================================================================

mesh <- function(x, y, z = NULL) { # x, y, z are vectors
  if (is.null(z)) 
    Mgrid.matrix(x, y)
  else  
    Mgrid.array(x, y, z)
}

Mgrid.matrix <- function(x, y) { 
   Nx <- length(x)
   Ny <- length(y)
   list(
    x = matrix(nrow = Nx, ncol = Ny, data = x),
    y = matrix(nrow = Nx, ncol = Ny, data = y, byrow = TRUE)
   )
}

Mgrid.array <- function(x, y, z) { 
   Nx <- length(x)
   Ny <- length(y)
   Nz <- length(z)
   list(
    x = array(dim = c(Nx, Ny, Nz), data = x),
    y = aperm(array(dim = c(Ny, Nx, Nz), data = y), c( 2, 1, 3)),
    z = aperm(array(dim = c(Nz, Ny, Nx), data = z), c( 3, 2, 1))
   )
}


## =============================================================================
## images in 3-D
## =============================================================================

createsurfs <- function (x, y, z, colvar, names = c("x", "y", "z"), resfac) {  

  if (is.null (x) & is.null(colvar) |
      is.null (y) & is.null(colvar))
    stop ("'colvar' cannot be NULL if ", names[1], " and/or ", names[2], " are NULL")
  
  if (is.null (x))
    x <- seq(0, 1, length.out = nrow(colvar))

  if (is.array(x)) {
    if (length(dim(x)) == 1)
      x <- as.vector(x)
    else if (length(dim(x)) == 2)
      x <- as.matrix(x)
  }
  
  if (! is.null(colvar)) {
    if (is.vector(x) & length(x) != nrow(colvar))
      stop (names[1], " should be a vector of length = nrow(colvar) or be NULL")
    else if (is.matrix(x))
      if (any(dim(x) - dim(colvar) < 0) | any(dim(x) - dim(colvar) > 1))
        stop(names[1], " not compatible with 'colvar'")
  }
      
  if (! is.matrix(x) & ! is.vector(x))
    stop(names[1], " should be a vector, a matrix or one value")
              
  if (is.null (y))
    y <- seq(0, 1, length.out = ncol(colvar))
  
  if (is.array(y)) {
    if (length(dim(y)) == 1)
      y <- as.vector(y)
    else if (length(dim(y)) == 2)
      y <- as.matrix(y)
  }
  if (! is.null(colvar)) {
    if (is.vector(y) & length(y) != ncol(colvar))
      stop (names[2], " should be a vector of length = ncol(colvar) or be NULL")
    else if (is.matrix(y))
      if (any(dim(y) - dim(colvar) < 0) | any(dim(y) - dim(colvar) > 1))
        stop(names[2], "not compatible with 'colvar'")
  }
  
  if (! is.matrix(y) & ! is.vector(y))
    stop(names[2], "should be a vector, a matrix or one value")

  if (is.vector(x) & ! is.vector(y) | 
      is.vector(y) & ! is.vector(x) )
    stop (names[1], " and ", names[2], " should both be a vector or both a matrix")

  if (any(resfac != 1)) {   
    res <- changeres(resfac, x, y, colvar)
    x <- res$x
    y <- res$y
    colvar <- res$z
  }

  if (is.vector(x)) {
    M <- mesh(x, y)
    x <- M$x
    y <- M$y
  }  
  z <- matrix(nrow = nrow(x), ncol = ncol(x), data = z)  
  list(A = x, B = y, C = z, colvar = colvar)
}

## =============================================================================
## main function
## =============================================================================

image3D <- function(x = NULL, y = NULL, z = NULL, ..., 
                  colvar = NULL, phi = 40, theta = 40,
                  col = NULL,  NAcol = "white", breaks = NULL,
                  border = NA, facets = TRUE,
                  colkey = NULL, resfac = 1,
                  panel.first = NULL,
                  clim = NULL, clab = NULL, bty = "b",
                  inttype = 1, add = FALSE, plot = TRUE){

  isconstant <- NULL
  if (length(x) == 1)
    isconstant <- c(isconstant, 1)
  if (length(y) == 1)
    isconstant <- c(isconstant, 2)
  if (length(z) == 1)
    isconstant <- c(isconstant, 3)
  if (length(isconstant) != 1)
    stop ("one of the values 'x' 'y', or 'z' should be one value")

  if (isconstant == 3) {
    ll <- createsurfs(x, y, z, colvar, c("x", "y", "z"), resfac)
    x <- ll$A
    y <- ll$B
    z <- ll$C 
  } else if (isconstant == 1) {
    ll <- createsurfs(y, z, x, colvar, c("y", "z", "x"), resfac)
    y <- ll$A
    z <- ll$B
    x <- ll$C 
  } else {
    ll <- createsurfs(x, z, y, colvar, c("x", "z", "y"), resfac)
    x <- ll$A
    z <- ll$B
    y <- ll$C 
  }
  if (any(resfac != 1))    
    colvar <- ll$colvar

  pmat <- surf3D (x = x, y = y, z = z, colvar = colvar, 
                  phi = phi, theta = theta,
                  col = col,  NAcol = NAcol, breaks = breaks,
                  border = border, facets = facets,
                  colkey = colkey, panel.first = panel.first,
                  clim = clim, clab = clab, bty = bty,
                  inttype = inttype, add = add, plot = plot, ...)
  
  invisible (pmat)
}                  


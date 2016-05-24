## =============================================================================
## 3-D visualisation of volumetric data using points
## =============================================================================
# x, y, z vectors, colvar: array

voxel3D <- function(x, y, z, colvar, ..., 
                    phi = 40, theta = 40, 
                    level = mean(colvar, na.rm = TRUE), eps = 0.01,
                    operator = "=", col = NULL,
                    NAcol = "white", breaks = NULL, colkey = FALSE,
                    panel.first = NULL, bty = "b", 
                    add = FALSE, plot = TRUE) {
  plist <- initplist(add)

  dot <- splitdotpersp(list(...), bty, NULL, x, y, z, plist = plist, breaks = breaks)

  if (length(level) != 1 & operator != "<>" )
    stop ("'level' should be one number if 'operator' not equal to '<>'")
  else if (length(level) != 2 & operator == "<>" )
    stop ("'level' should be two numbers if 'operator' equals '<>'")
  
  if (is.null(plist)) {
    do.call("perspbox", c(alist(x = range(x), y = range(y), 
             z = range(z, na.rm = TRUE),
             phi = phi, theta = theta, colkey = colkey, plot = plot,
             col = col), dot$persp))
    plist <- getplist()
  }  
  if (is.function(panel.first)) 
    panel.first(plist$mat)         

  vox <- createvoxel (x, y, z, colvar, level, eps, operator)
  if (is.null(vox)) {
    rc <- range(colvar, na.rm = TRUE)
    stop("no points selected - change 'level': valid range ", formatC(rc[1]), 
      " - ", formatC(rc[2]))
  }
  if (operator == "=")
    colvar <- vox$z
  else 
    colvar <- vox$cv    
  
  do.call("scatter3D", c(alist(x = vox$x, y = vox$y, z = vox$z, 
          add = TRUE, col = col, NAcol = NAcol, breaks = breaks,
          colkey = FALSE, plot = plot, alpha = dot$alpha), dot$points))
  plist <- getplist()
  invisible(plist$mat)
}

## =============================================================================

createvoxel <- function (x, y, z, colvar, level = mean(colvar, na.rm = TRUE), 
                         eps = 0.01, operator = "=") {

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

  clim <- range(colvar, na.rm = TRUE)     
  eps <- diff(clim) * eps
  
  if (operator == "=")
    ijk <- which(abs(colvar - level) < eps , arr.ind = TRUE)
  else if (operator == "<")
    ijk <- which(colvar < level, arr.ind = TRUE)
  else if (operator == ">")
    ijk <- which(colvar > level, arr.ind = TRUE)
  else if (operator == "<>")
    ijk <- which(colvar > level[1] & colvar < level[2], arr.ind = TRUE)
  else
    stop ("'operator' should be one of '=', '<', or '>'")

  if (nrow(ijk) > 1) 
    list (x = x[ijk[,1]], y = y[ijk[,2]], z = z[ijk[,3]], cv = colvar[ijk])
}

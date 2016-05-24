## =============================================================================
## Contours in 3-D
## =============================================================================

createsegms <- function (x, y, z, colvar, names = c("x", "y", "z"), 
  dot, col, clim, dDepth, plist, levels, addbox) {  

  contour <- list(args = dot$points) 

  contour$args$col <- col
  contour$args$levels <- levels

  if (!ispresent(colvar)) 
    stop ("'colvar' should be present for contour3D")

  if (is.null(x))
    x <- seq(0, 1, length.out = nrow(colvar))
  if (is.null(y))
    y <- seq(0, 1, length.out = ncol(colvar))

  if (is.array(x) & length(dim(x)) == 1)
    x <- as.vector(x)
  if (is.array(y) & length(dim(y)) == 1)
    y <- as.vector(y)

  if (! is.vector(x))
    stop(names[1], " should be a vector")
    
  if (length(x) != nrow(colvar))
    stop (names[1], " should be a vector of length = nrow(colvar) or be NULL")

  if (! is.vector(y))
    stop(names[2], " should be a vector")

  if (length(y) != ncol(colvar))
    stop (names[2], " should be a vector of length = ncol(colvar) or be NULL")

  if (! is.matrix(z)) {
    if (length(z) > 1)
      stop("'z'  should be a matrix or one value") 
    contour$side <- z
    z <- colvar 
  } else {
    if (length(x) != nrow(z))
      stop(names[1], " should be of length = nrow(",names[3],")")
    if (length(y) != ncol(z))
      stop(names[2], " should be of length = ncol(",names[3],")")
    contour$side <- "z"
    
  }
  
 # create contours
  segm <- contourfunc(contour, x, y, z, plist, cv = colvar, 
    clim = clim, dDepth = dDepth, addbox = addbox) 
  names.from <- paste(names, ".from", sep = "")
  names.to <- paste(names, ".to", sep = "")

  names(segm)[1:6] <- c(names.from[1], names.to[1], 
                        names.from[2], names.to[2], 
                        names.from[3], names.to[3])
  return(segm)
}

## =============================================================================
## main function
## =============================================================================

contour3D <- function(x = NULL, y = NULL, z = NULL, ..., 
                  colvar = NULL, phi = 40, theta = 40,
                  col = NULL, colkey = NULL, 
                  panel.first = NULL,
                  clim = NULL, clab = NULL, bty = "b",
                  dDepth = 1e-1, addbox = TRUE,
                  add = FALSE, plot = TRUE){

  xlim <- ylim <- zlim <- c(0, 1)
  if (!is.null(x))
    xlim <- range(x, na.rm = TRUE)
  if (!is.null(y))
    ylim <- range(y, na.rm = TRUE)
  if (!is.null(z))
    zlim <- range(z, na.rm = TRUE)

  plist <- initplist(add)

  dot <- splitdotpersp(list(...), bty, FALSE, xlim, ylim, zlim,
    plist = plist, breaks = NULL)

  dots <- dot$points    
  levels <- dots$levels
  dots$levels <- NULL
  
  if (is.null(clim)) 
    clim <- range(colvar, na.rm = TRUE)

  clog <- FALSE
  if (! is.null(dots$log)) {
    if (length(grep("c", dots[["log"]])) > 0) {
      dots[["log"]] <- gsub("c", "", dots[["log"]])
      clog <- TRUE
    }
    if (dots[["log"]] == "")
      dots[["log"]] <- NULL
  }
  if (clog) {                
    colvar <- log(colvar)
    clim <- log(clim)
  }

  if (is.null(levels)) {
    nlevs <- dots$nlevels

    if (is.null(nlevs))
      nlevs <- 10

    if (clog) 
      levels <- exp(pretty(log(clim), nlevs))
    else
      levels <- pretty(clim, nlevs)
  }
  nlevs <- length(levels)

  if (is.null(col))
    col <- jet.col(nlevs)
  if (! is.null(dot$alpha)) 
    col <- setalpha(col, dot$alpha)

  iscolkey <- is.colkey(colkey, col)
  if (iscolkey) 
    colkey <- check.colkey(colkey)
   
  if (is.null(plist)) {
    do.call("perspbox", c(alist(xlim, ylim, zlim,  
                     phi = phi, theta = theta, plot = plot, 
                     colkey = colkey, col = col), dot$persp))
    plist <- getplist()
  }  
  if (is.function(panel.first)) 
    panel.first(plist$mat)         
                                 
  isconstant <- NULL
  ismatrix <- NULL
  if (length(x) == 1)
    isconstant <- c(isconstant, 1)
  else if (is.matrix(x))
    ismatrix <- c(ismatrix, 1)
  
  if (length(y) == 1)
    isconstant <- c(isconstant, 2)
  else if (is.matrix(y))
    ismatrix <- c(ismatrix, 2)

  if (length(z) == 1)
    isconstant <- c(isconstant, 3)
  else if (is.matrix(z))
    ismatrix <- c(ismatrix, 3)

  if (length(isconstant) > 1)
    stop ("only one of the values 'x' 'y', or 'z' can be one value")

  if (length(ismatrix) > 1)
    stop ("only one of the values 'x' 'y', or 'z' can be a matrix")

  if (length(isconstant) > 1)
    stop ("only one of the values 'x' 'y', or 'z' can be one value")

  if (length(ismatrix) == 0 & length(isconstant) == 0)
    stop ("exactly one of the values 'x' 'y', or 'z' should be a matrix or one value")

  ismapped <- c(isconstant, ismatrix)
  if (ismapped == 3) {
    segm <- createsegms(x, y, z, colvar, c("x", "y", "z"), dot, 
      col, clim, dDepth, plist, levels, addbox)
  } else if (ismapped == 1) {
    segm <- createsegms(y, z, x, colvar, c("y", "z", "x"), dot, 
      col, clim, dDepth, plist, levels, addbox)
  } else {
    segm <- createsegms(x, z, y, colvar, c("x", "z", "y"), dot, 
      col, clim, dDepth, plist, levels, addbox)
  }

  if (iscolkey) {
    if (is.null(colkey$at))
      colkey$at <- levels

    if (nlevs > 1) {
     # for colors: 
      dz <- c(-diff(levels[1:2]), diff(levels[(nlevs-1):nlevs])) * 0.5
      clim <- range(levels) + dz
    }
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
      dot$clog, type = "contour3D", breaks = NULL)
  }
  plist <- plot.struct.3D(plist, segm = segm, plot = plot)  

  setplist(plist)   
  invisible(plist$mat)
  
}                  


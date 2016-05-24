moviepersp <- function(z, x = NULL, y = NULL, t = NULL,
    colvar = z, tdim = 1,         # 1st dimension of z is time
    col = jet.col(100), NAcol = "white", breaks = NULL,
    colkey = FALSE, clim = NULL, clab = NULL, wait  = NULL, ask = FALSE,
     add = FALSE, basename = NULL, ... ) {

  dZ <- dim(z)
  if (length(dZ) != 3)
    stop("'z' should be a 3-dimensional array")
  if (tdim > 3 | tdim < 1)
    stop ("'tdim' should refer to the position in 'z' of 't', either 1, 2, 3")
  dXY <- dZ[-tdim]
  if (is.null(x))
    x <- seq(0, 1, length.out = dXY[1])
  else {
   if (! is.vector(x))
     stop("'x' should be a vector with length = dim(z)[-tdim][1]")

   if (length(x) != dXY[1])
     stop("'x' should be a vector with length = dim(z)[-tdim][1]")
  }

  if (is.null(y))
    y <- seq(0, 1, length.out = dXY[2])
  else {
   if (! is.vector(y))
     stop("'y' should be a vector with length = dim(z)[-tdim][2]")

   if (length(y) != dXY[2])
     stop("'y' should be a vector with length = dim(z)[-tdim][2]")
  }

  if (is.null(t))
    t <- seq(0, 1, length.out = dZ[tdim])
  else {
   if (! is.vector(t))
     stop("'t' should be a vector with length = dim(z)[tdim]")

   if (length(t) != dZ[tdim])
     stop("'t' should be a vector with length = dim(z)[tdim]")
  }

  if (is.null(col) & is.null(breaks))
    col <- jet.col(100)
  else if (is.null(col))
    col <- jet.col(length(breaks)-1)
  breaks <- check.breaks(breaks, col)

  if (length(col) == 1)
    colvar <- NA

  if (is.null(colvar))
    colvar <- z
  else {
    dC <- dim(colvar)
    if (any(dC - dZ != 0))
      stop ("'colvar' should be of same dimension as 'z'")
  }
   if (is.null(clim))  {
    clim <- range(colvar, na.rm = TRUE)
  } else {
    colvar[colvar < min(clim)] <- NA
    colvar[colvar > max(clim)] <- NA
  }

  Col <- variablecol (colvar, col, NAcol, clim, breaks)
  dim(Col) <- dim(colvar)
  dots <- list(...)

  if (tdim == 1)
    getval <- function(A, i)  A[i,,]
  else if (tdim == 2)
    getval <- function(A, i)  A[,i,]
  else if (tdim == 3)
    getval <- function(A, i)  A[,,i]
  else
    stop ("'tdim' should be 1, 2, or 3")

  if (!interactive())
    ask  <- FALSE

  new <- dots$new
  dots$new <- NULL
  if (is.null(new))
    new <- TRUE

  fi <- 0

  if (is.null(dots$zlim))
    zlim <- range(z, na.rm = TRUE)
  else
    zlim <- dots$zlim
  dots$zlim <- NULL

  cv <- colvar
  Main <- dots$main
  dots$main <- NULL

   if (add)
    plist <- getplist()

  for (i in 1:length(t)) {
    time <- t[i]
    main <- paste(Main, time)
     if (is.array(colvar))
       cv <- getval(colvar, i)
     else
       cv <- colvar
    if (add) {
      setplist(plist)
      plotdev()
    }
    do.call("persp3D", c(alist(persp3D(x, y, z = getval(z, i), colvar = cv,
       add = add, col = col, breaks = breaks, clim = clim,
       main = main, zlim = zlim, colkey = colkey), dots)))
    if (ask)
     readline("Hit enter to change surface")
    else if (! is.null(wait))
     Sys.sleep(wait)
    if (! is.null(basename)) {
      filename <- paste(basename, formatC(as.integer(fi), width = 4, format = "d", flag = "0"),".png", sep="")
      rgl.snapshot(filename)
      fi <- fi + 1
    }
  }
}


## =============================================================================
## =============================================================================
## movies of 3D perspective plot 
## =============================================================================
## =============================================================================

moviepersp3D <- function(z, x = NULL, y = NULL, t = NULL, 
    colvar = z, tdim = 1,         # 1st dimension of z is time
    col = jet.col(100), NAcol = "white", breaks = NULL,
    colkey = FALSE, clim = NULL, clab = NULL, wait  = NULL, ask = FALSE,
     add = FALSE, basename = NULL, ... ) {

  dZ <- dim(z)
  if (length(dZ) != 3)
    stop("'z' should be a 3-dimensional array")
  if (tdim > 3 | tdim < 1)
    stop ("'tdim' should refer to the position in 'z' of 't', either 1, 2, 3")
  dXY <- dZ[-tdim]
  if (is.null(x))
    x <- seq(0, 1, length.out = dXY[1])
  else {
   if (! is.vector(x))
     stop("'x' should be a vector with length = dim(z)[-tdim][1]")

   if (length(x) != dXY[1])
     stop("'x' should be a vector with length = dim(z)[-tdim][1]")
  }
  
  if (is.null(y))
    y <- seq(0, 1, length.out = dXY[2])
  else {
   if (! is.vector(y))
     stop("'y' should be a vector with length = dim(z)[-tdim][2]")

   if (length(y) != dXY[2])
     stop("'y' should be a vector with length = dim(z)[-tdim][2]")
  }

  if (is.null(t))
    t <- seq(0, 1, length.out = dZ[tdim])
  else {
   if (! is.vector(t))
     stop("'t' should be a vector with length = dim(z)[tdim]")

   if (length(t) != dZ[tdim])
     stop("'t' should be a vector with length = dim(z)[tdim]")
  }

  if (is.null(col) & is.null(breaks))
    col <- jet.col(100)
  else if (is.null(col))
    col <- jet.col(length(breaks)-1)
  breaks <- check.breaks(breaks, col)

  if (length(col) == 1)
    colvar <- NA

  if (is.null(colvar))
    colvar <- z
  else {
    dC <- dim(colvar)
    if (any(dC - dZ != 0))
      stop ("'colvar' should be of same dimension as 'z'") 
  }  
   if (is.null(clim))  {
    clim <- range(colvar, na.rm = TRUE)
  } else {
    colvar[colvar < min(clim)] <- NA
    colvar[colvar > max(clim)] <- NA
  }

  Col <- variablecol (colvar, col, NAcol, clim, breaks)
  dim(Col) <- dim(colvar)
  dots <- list(...)

  if (tdim == 1)
    getval <- function(A, i)  A[i,,]
  else if (tdim == 2)
    getval <- function(A, i)  A[,i,]
  else if (tdim == 3)
    getval <- function(A, i)  A[,,i]
  else
    stop ("'tdim' should be 1, 2, or 3")

  popit <- FALSE
  if (!interactive()) 
    ask  <- FALSE
  
  new <- dots$new
  dots$new <- NULL
  if (is.null(new))
    new <- TRUE 

  fi <- 0
  
  if (is.null(dots$zlim))
    zlim <- range(z, na.rm = TRUE)
  else
    zlim <- dots$zlim
  dots$zlim <- NULL

  cv <- colvar
  Main <- dots$main
  dots$main <- NULL
  idsmain <- NULL
  if (!add) {
     do.call("persp3Drgl", c(alist(x = range(x), y = range(y), z = diag(zlim),
       col = col, breaks = breaks, add = add, clim = clim, new = new,
       main = "", zlim = zlim, colkey = colkey), dots))
       rglids <- rgl.ids()
       rgl.pop(id = rglids[1,1])
#     on.exit(par3d(save))
  idsmain <- getplist()$rgl$D$main
  }
  for (i in 1:length(t)) {
    time <- t[i]
    main <- paste(Main, time)
     if (is.array(colvar))
       cv <- getval(Col, i)
    if (popit)
      rgl.pop(id = 0)

    if (! is.null(idsmain))
        rgl.pop(type = "shapes", id = idsmain)
    idsmain <- mtext3d(main, "x++", line = 4)
    persp3d(x, y, z = getval(z, i), col = cv, add = TRUE, #aspect = FALSE,
            front = "filled", back = "filled")#, smooth = smooth, alpha = fb$alpha)

    popit <- TRUE
    if (i < length(t)) {
    if (ask)
     readline("Hit enter to change surface")
    else if (! is.null(wait))
     Sys.sleep(wait)
    }
    if (! is.null(basename)) {
      filename <- paste(basename, formatC(as.integer(fi), width = 4, format = "d", flag = "0"),".png", sep="")
      rgl.snapshot(filename) 
      fi <- fi + 1
    }
  }
}


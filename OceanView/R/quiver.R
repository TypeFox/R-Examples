
## =============================================================================
## =============================================================================
## QUIVER FUNCTIONS
## =============================================================================
## =============================================================================

checkinput <- function(u, v, x = NULL, y = NULL, scale = 1, by = 1, 
  xlim = NULL, ylim = NULL, maxspeed = NULL, Log = FALSE, plist = NULL, 
  add = FALSE) {  
  
  if (is.null(x)) 
    x <- seq(0, 1, length.out = nrow(u))
  if (is.null(y)) 
    y <- seq(0, 1, length.out = ncol(v))

  if (is.array(x)) {
    if (length(dim(x)) == 1) 
      x <- matrix(ncol = ncol(u), nrow = length(x), data = x)
    else if (length(dim(x)) > 2) 
      stop ("'x' cannot be an array with more than 2 dimension")
  }

  if (is.array(y)) {
    if (length(dim(y)) ==1) 
      y <- matrix(nrow = nrow(v), ncol = length(y), data = y, byrow = TRUE)
    else if (length(dim(y)) > 2) 
      stop ("'y' cannot be an array with more than 2 dimension")
  }

  if (is.vector(x)) 
    x <- matrix(ncol = ncol(u), nrow = length(x), data = x)
  if (is.vector(y)) 
    y <- matrix(nrow = nrow(v), ncol = length(y), data = y, byrow = TRUE)
  
  div <- dim(u) - dim(v)
  
# u and v may have a dimension 1 different, if they are defined on the edges.       
  if (div[1] == -1) 
    v <- 0.5*(v[-1,] + v[-nrow(v),])
  else if (div[1] == 1) 
    u <- 0.5*(u[-1,] + u[-nrow(u),])
  else if (div[1] != 0)
    stop("dimensions of 'u' and 'v' not compatible") 
    
  if (div[2] == 1) 
    u <- 0.5*(u[,-1] + u[,-ncol(u)])
  else if (div[2] == -1)  
    v <- 0.5*(v[,-1] + v[,-ncol(v)])
  else if (div[2] != 0)
    stop("dimensions of 'u' and 'v' not compatible") 

# Check x- and y
  Dim <- dim(u)
    
  if (is.matrix(x)) {
    if (ncol(x) == Dim[2] +1)
      x <- 0.5*(x[,-1] + x[,-ncol(x)])
    if (nrow(x) == Dim[1] +1)
      x <- 0.5*(x[-1,] + x[-nrow(x),])
    if (ncol(x) != Dim[2] | nrow(x) != Dim[1])
      stop ("'x' not compatible with u or v")
  }
  if (is.matrix(y)) {
    if (ncol(y) == Dim[2] +1)
      y <- 0.5*(y[,-1] + y[,-ncol(y)])
    if (nrow(y) == Dim[1] +1)
      y <- 0.5*(y[-1,] + y[-nrow(y),])
    if (ncol(y) != Dim[2] | nrow(y) != Dim[1])
      stop ("'y' not compatible with u or v")
  }

  if (is.null(x)) 
    x <- (row(u)-0.5)/nrow(u)
  if (is.null(y)) 
    y <- (col(v)-0.5)/ncol(v)

# ------------------------------------------------------------------------------
# select the elements of x, y, u and v
# ------------------------------------------------------------------------------    
  if (! is.null(by)) {                     
    by <- rep(by, length = 2)
    ix <- seq(1, nrow(u), by = by[1])
    iy <- seq(1, ncol(u), by = by[2])
    u <- u[ix,iy]
    x <- x[ix,iy]
    v <- v[ix,iy]
    y <- y[ix,iy]
  } else {
    ix <- 1:nrow(u)
    iy <- 1:ncol(u)
  }
  if (! is.null(xlim)) {
    ii <- which (x < min(xlim) | x > max(xlim))
    u[ii] <- v[ii] <- 0
  }
  if (! is.null(ylim)) {
    ii <- which (y < min(ylim) | y > max(ylim))
    u[ii] <- v[ii] <- 0
  }
  isna <- which (is.na(u))
  u[is.na(u)] <- 0
  v[is.na(v)] <- 0
    
# ------------------------------------------------------------------------------
# size of the arrows
# ------------------------------------------------------------------------------    
  speed <- sqrt(u^2 + v^2)
  
  if (is.null(maxspeed))
    maxspeed <- max(speed) 
  else {
    if (maxspeed <= 0)
      stop ("'speed.max' should be >= 0")
    speed <- pmin(speed, maxspeed)
  }  
  if (maxspeed == 0) 
    maxspeed <- 1e-16
  if (!is.null(plist$quiver)) {
    xr <- plist$quiver$xr
    yr <- plist$quiver$yr
    maxspeed <- plist$quiver$maxspeed
    if (Log) maxspeed <- exp(maxspeed)
  } else if (add) {
    pusr <- par("usr")
    xr <- diff (pusr[1:2]) / max(dim(u))
    yr <- diff (pusr[3:4]) / max(dim(u))
  
  } else {
  xr <- diff (range(x)) / max(dim(u))
  yr <- diff (range(y)) / max(dim(u))
  }
  if (!is.null(scale)) {
    u <- u * scale / maxspeed * xr 
    v <- v * scale / maxspeed * yr 
  }
  if (Log) speed <- log(speed)
  if (Log) maxspeed <- log(maxspeed)
  list(x = x, y = y, u = u, v = v, xr = xr, yr = yr,
       speed = speed, maxspeed = maxspeed, 
       isna = isna, ix = ix, iy = iy)
}

## =============================================================================
## Actual quiver functions
## =============================================================================

quiver2D <- function(u, ...) UseMethod ("quiver2D")
quiver2D.default <- function (u, ...) quiver2D.matrix(u, ...)

quiver2D.matrix  <- function(u, v, x = NULL, y = NULL, colvar = NULL, ..., 
                    scale = 1, arr.max = 0.2, arr.min = 0, speed.max = NULL,
                    by = NULL, type = "triangle", 
                    col = NULL, NAcol = "white", breaks = NULL, colkey = NULL,
                    mask = NULL, image = FALSE, contour = FALSE, 
                    clim = NULL, clab = NULL, add = FALSE, plot = TRUE)  {
# ------------------------------------------------------------------------------
# check input
# ------------------------------------------------------------------------------
  if (! is.null(mask) & add == TRUE)
    stop ("cannot combine a 'mask' with 'add  = TRUE'")

  dots  <- splitpardots( list(...) )
  dp    <- dots$points
  dm    <- dots$main 
 
  if (add) 
      plist <- getplist()
  else plist <- NULL
  setplist(plist)
 
  # colors and color variable
  if (! is.null(colvar)) {
    varlim <- clim
    if (is.null(varlim)) 
      varlim <- range(colvar, na.rm = TRUE)
  
    if (any (dim(colvar) - c(nrow(u), ncol(v)) != 0)) 
      stop ("dimension of 'colvar' not compatible with dimension of 'u' and 'v'")

    if (! is.null(by)) {                     
      by <- rep(by, length = 2)
      ix <- seq(1, nrow(u), by = by[1])
      iy <- seq(1, ncol(u), by = by[2])
      colvar <- colvar[ix, iy]
    }  

    if (is.null(col))
      col <- jet.col(100)

    if (dots$clog) {
      colvar <- log(colvar)
      if (! is.null(clim)) 
        clim <- log(clim) 
    }
    
    iscolkey <- is.colkey(colkey, col)    # check if colkey is needed
    if (iscolkey) {
      colkey <- check.colkey(colkey)
      if (! add) 
        plist$plt$main <- colkey$parplt
    }  
    par (plt = plist$plt$main)
    
    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
    
    Col <- variablecol(colvar, col, NAcol, clim, breaks) # generate color scheme

    pltori <- plist$plt$ori
  } else {
    Col <- col
    if (is.null(Col)) 
      Col <- "black"
    iscolkey <- FALSE
  }   
  if (plot) {
    image   <- check.args(image, NULL)
    contour <- check.args(contour, NULL)

  # images, contours or masked cells
    if (! is.null(mask)) {
      if (image$add)
        stop ("cannot have both 'image' and 'mask' specified")

      maskNAcol <- NULL
      X <- x; Y <- y
      if (is.list(mask)) {
        maskNAcol <- mask$NAcol  
        if (! is.null(mask$x)) 
          X <- mask$x 
        if (! is.null(mask$y)) 
          Y <- mask$y 
        if (! is.null(mask$z)) 
          mask <- mask$z 
      }
      if (is.null(maskNAcol))
        maskNAcol <- "black"
      mask[!is.na(mask)] <- 0
      mask[is.na(mask)]  <- 1
      do.call ("image2D", c(alist(z = mask, x = X, y = Y, colkey = FALSE,
                 col = c("white", maskNAcol)), dm))
      add <- TRUE
    }

    if (image$add) {
      if (is.null(image$args$z))
        stop("image$'z' should be specified if image is a list")
      if (is.null(image$args$col))
        image$args$col <- jet.col(100)
      if (!add) 
        image$args <- c(image$args, dm)
      do.call("image2D", c(alist(colkey = image$colkey), image$args))   
      add <- TRUE
    }

    if (contour$add) {
     if (add) 
      contour$args$colkey <- FALSE
      if (is.null(contour$args$z))
        stop("contour$'z' should be specified if 'contour' is a list")

      if (!add) 
        contour$args <- c(contour$args, dm)

      do.call("contour2D", c(alist(add = add), contour$args))
      add <- TRUE
    }
  } # plot
  # karline: log of arrows log = "a""
  Log <- FALSE
  if (! is.null(dm$log)) {
    if (length(grep("a", dm[["log"]])) > 0) {
      dm[["log"]] <- gsub("a", "", dm[["log"]])
      Log <- TRUE
      if (dm[["log"]] == "")
        dm[["log"]] <- NULL
    }
  }

  MM <- checkinput(u, v, x, y, scale, by = by, xlim = dm$xlim, 
    ylim = dm$ylim, maxspeed = speed.max, Log, plist, add) 
  x <- MM$x
  y <- MM$y
  u <- MM$u
  v <- MM$v
   
  xto <- x + u
  yto <- y + v

 # transpose dp elements that are matrices
  
  dp <- lapply(dp, FUN = function(x) 
                    if (is.matrix(x)) x[MM$ix, MM$iy] else x)
  
#  dp$arr.length <- MM$speed / MM$maxspeed * (arr.max - arr.min) + arr.min
  dp$length <- MM$speed / MM$maxspeed * (arr.max - arr.min) + arr.min

  if (plot) {
    if (! add) {
      if (is.null(dm$xlab)) 
        dm$xlab <- "x"
      if (is.null(dm$ylab)) 
        dm$ylab <- "y"
    }
  
    dp$arr.type <- NULL 
    if (is.null(dp$lwd))   
      dp$lwd <- 1
  
    do.call("arrows2D", c(alist(x, y, xto, yto, col = Col, type = type, add = add),  dm, dp))
  
    if (iscolkey) {
      colkey$parleg <- colkey$parplt <- NULL    
      do.call("colkey", c(alist(col = col, clim = varlim, clab = clab, 
        clog = dots$clog, add = TRUE), colkey))
      par(plt = pltori)  
    }    
    par(mar = par("mar")) 

  } #plot
  if (Log) MM$maxspeed <- exp(MM$maxspeed)
  if (! add) {
    plist <- getplist()
    plist$quiver <- list(xr = MM$xr, yr = MM$yr, maxspeed = MM$maxspeed)
    setplist(plist)
  }
    
  invisible(list(x0 = x, y0 = y, x1 = xto, y1 = yto, col = Col, length = dp$length, 
    speed.max = MM$maxspeed))
}

## =============================================================================

quiver2D.array <- function (u, v, margin = c(1, 2), subset, ask = NULL, ...) {

  DD <- dim(u)
  if (length(DD) != 3)
    stop ("Can only make quiver of 3-D array, 'u' has dimension ", length(DD))
  if (length(dim(v)) != 3)
    stop ("Can only make quiver of 3-D array, 'v' has dimension ", length(dim(v)))

  if (length(margin) != 2)
    stop ("'margin' should contain two numbers, the x, y subscripts of which to make quiver plots")
   
  if ( max(margin) > 3 | min (margin) < 1)
    stop ("indexes in 'margin' should be inbetween 1 and 3")

  index <- (1:3) [- margin]

  if (index > 3 || index <1)
    stop ("'index' to loop over should be inbetween 1 and 3")
  
  x <- 1:DD[index]
  ldots  <- list(...)
  colvar <- ldots$colvar 

  if (!missing(subset)){
     if (is.numeric(subset)) { 
       isub <- subset
     } else {
      e <- substitute(subset)
      r <- eval(e, as.data.frame(x), parent.frame())
      if (!is.logical(r))
          stop("'subset' must evaluate to logical")
      isub <- r & !is.na(r)
      isub <- which(isub)
      if (length(isub) == 0)
        stop("cannot continue: nothing selected - check 'subset'")
     }   
    
  } else isub <- x

  np     <- length(isub)
 # Set par mfrow and ask
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (is.null(ldots$main)) 
    title <- isub
  else   
    title <- rep(ldots$main, length.out = length(isub))
  i1 <- 1
  for (i in isub) {
    LL <- ldots
    if (index == 1) { 
      LL$u <- u[i, , ] 
      LL$v <- v[i, , ] 
      if (! is.null(colvar))  
        LL$colvar <- colvar[i, , ] 
    } else if (index == 2) {
      LL$u <- u[ ,i , ] 
      LL$v <- v[ ,i , ] 
      if (! is.null(colvar))  
        LL$colvar <- colvar[, i, ] 
    } else {
      LL$u <- u[ ,, i ]
      LL$v <- v[ ,, i ]
      if (! is.null(colvar))  
        LL$colvar <- colvar[, , i] 
    }  
    if (margin[2] < margin[1]) {
      LL$u <- t(LL$u)
      LL$v <- t(LL$v)
    }  
    LL$main <- title[i1]
    i1 <- i1 +1
    do.call(quiver2D, LL)
    #quiver(u = uu, v = vv, ...) 
  }
  
}

## =============================================================================
## =============================================================================
## flowpaths based on flow fields
## =============================================================================
## =============================================================================

flowpath <- function (u, v, x = NULL, y = NULL, startx = NULL, starty = NULL, ...,
                       scale = 1,  numarr = 0, arr.length = 0.2, 
                       maxstep = 1000, add = FALSE, plot = TRUE) {

  if (is.null(x)) 
    x <- seq(0, 1, length.out = nrow(u))
  if (is.null(y)) 
    y <- seq(0, 1, length.out = ncol(v))

  MM <- checkinput(u, v, x, y, scale, by = 1) 
  
  x <- MM$x
  y <- MM$y
  u <- MM$u
  v <- MM$v

  Nx <- nrow(x)
  Ny <- ncol(y)

  xx <- x[,1]
  yy <- y[1,]
  dx <- c(diff(xx), xx[Nx]-xx[Nx-1])
  dy <- c(diff(yy), yy[Ny]-yy[Ny-1])

  if (is.null(startx) ) {
    startx <- c(rep(xx[1], Ny),rep(xx[Nx], Ny), xx, xx) 
    starty <- c(yy, yy, rep(yy[1], Nx), rep(yy[Ny], Nx))
  } 
  
  lx <- length(startx)
  ly <- length(starty)
  if (lx != ly) {
    startx <- rep(startx, len = ly)
    starty <- rep(starty, len = length(startx))
  }
  
  xto <- startx
  yto <- starty
  
  if (min(xto) < min(x)) 
    stop("'x' should embrace 'startx'")
  if (min(yto) < min(y)) 
    stop("'y' should embrace 'starty'")
    
  interp <- function (xto, yto, u) {

  # find embracing values : first interval
    ix <- FindInterval(xto, xx)
    iy <- FindInterval(yto, yy)

  # next interval
    ixp1 <- pmin(ix+1,Nx)
    iyp1 <- pmin(iy+1,Ny)

  # interpolation factor
    xfac <- (xto-xx[ix])/dx[ix]
    yfac <- (yto-yy[iy])/dy[iy]

  # interpolate
    (1-yfac)*((1-xfac)*u[cbind(ix,iy)]+xfac*u[cbind(ixp1,iy)]) +
    yfac*((1-xfac)*u[cbind(ix,iyp1)]+xfac*u[cbind(ixp1,iyp1)])

  } # end Transf

  if (plot) {
    dotlist  <- splitpardots( list(...) )
    dp <- dotlist$points
    dm <-  dotlist$main

    if (is.null(dm$xlab)) 
      dm$xlab <- "x"
    if (is.null(dm$ylab)) 
      dm$ylab <- "y"
 
    if (!add) 
      do.call("matplot", c(alist(x, y), type = "n", dm))
    dpa <- dp                   
    dp[grep("arr",names(dp))] <- NULL
    if ( numarr >0) {
      if (is.null(dpa$arr.type))  dpa$arr.type <- "triangle"  
      dpa$arr.length <- arr.length
    }
  }
#  XY <- cbind 
  xr <- range(xx)
  yr <- range(yy)
  mm <- min(diff(xx), diff(yy))*1e-12
  ij <- length(startx)

  xy <- XY <- NULL
  
  for (j in 1:ij) {
    xto <- startx[j]
    yto <- starty[j]       
    xy <- rbind(xy, c(xto, yto))
    for (i in 2: maxstep) {
      u2 <- interp(xto, yto, u) 
      v2 <- interp(xto, yto, v) 
      xto <- xto + u2
      yto <- yto + v2
      xy <- rbind(xy, c(xto, yto))

      if (xto < xr[1] | xto > xr[2] | yto < yr[1] | yto > yr[2] | 
          (abs(u2) < mm & abs(v2) < mm)) break
    }
    if (plot) {
    
      do.call ("lines", c(alist(xy[1:i, ]), dp))                     

      if (numarr > 0 & i > 2) {
        i1 <- 1/( numarr +1)
        ii <- pmax(2,as.integer(seq(i*i1, i *(1-i1), len =  numarr )))
        do.call("Arrows", 
                c(alist(xy[ii-1,1], xy[ii-1,2], xy[ii,1], xy[ii,2]), dpa, cex = 0.5))
      }
      XY <- rbind(XY, xy, c(NA, NA))
      xy <- NULL

    } else {
      XY <- rbind(XY, xy, c(NA, NA))
    }
  }
  invisible(XY)
}

## =============================================================================
## QUIVER using rgl graphics
## =============================================================================

quiver2Drgl <- function(u, v, x = NULL, y = NULL, colvar = NULL, ..., 
                    scale = 1, arr.max = 0.2, arr.min = 0, speed.max = NULL, 
                    by = NULL, type = "triangle", 
                    col = NULL, NAcol = "white", breaks = NULL,
                    mask = NULL, image = FALSE, contour = FALSE, 
                    colkey = NULL, clim = NULL, clab = NULL,
                    add = FALSE, plot = TRUE) {
  if (is.null(x))
    x <- seq(0, 1, length.out = nrow(u))
  xlim <- range(x)
  if (is.null(y))
    y <- seq(0, 1, length.out = ncol(v))
  ylim <- range(y)
  
  F <- quiver2D(u, v, x, y, colvar, scale = scale, 
              arr.max = arr.max, arr.min = arr.min, speed.max = speed.max, 
              by = by, plot = FALSE, col = col, NAcol = NAcol, breaks = breaks,
              colkey = colkey, clim = clim, clab = clab)

  arrows2Drgl(F$x0, F$y0, F$x1, F$y1, colvar = NULL, type = type, 
    col = F$col, NAcol = NAcol, add = add, code = 2, 
    length = F$length/2, angle = 15, ...)

  image   <- check.args(image, NULL)
  contour <- check.args(contour, NULL)
  addimage <- FALSE
 # images, contours or masked cells
  if (! is.null(mask)) {
    if (image$add)
      stop ("cannot have both 'image' and 'mask' specified")

      maskNAcol <- NULL
      X <- x; Y <- y
      if (is.list(mask)) {
        maskNAcol <- mask$NAcol  
        if (! is.null(mask$x)) 
          X <- mask$x 
        if (! is.null(mask$y)) 
          Y <- mask$y 
        if (! is.null(mask$z)) 
          mask <- mask$z 
      }
      mask[!is.na(mask)] <- 0
      mask[is.na(mask)]  <- 1
      if (is.null(maskNAcol)) 
        maskNAcol <- "black" 
      do.call ("image2Drgl", alist(z = mask, x = X, y = Y, dz = -0.01, colkey = FALSE,
                 plot = FALSE, add = TRUE, col = c("white", maskNAcol)))

      addimage <- TRUE
    }

    if (image$add) {
      image$args$colvar <- NULL 
      image$args$add <- image$args$plot <- NULL
      if (is.null(image$args$col))
        image$args$col <- jet.col(100)
      do.call("image2Drgl", c(alist(colkey = FALSE, plot = FALSE, add = TRUE, dz = -0.01),
          image$args))   
      addimage <- TRUE
    }

    if (contour$add) {
      contour$args$add <- contour$args$plot <- NULL
      contour$args$colvar <- NULL #contour$args$z
      do.call("contour2Drgl", c(alist(plot = FALSE, add = TRUE), contour$args, dz = 0.1))
      addimage <- TRUE
    }
}

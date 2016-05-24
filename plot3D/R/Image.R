
## =============================================================================
## =============================================================================
## image2D S3 functions   - this code can be improved.
## =============================================================================
## =============================================================================

image2D <- function(z, ...) UseMethod ("image2D")
image2D.default <- function (z, ...) image2D.matrix(z, ...)

## =============================================================================
## image2D function, input is a matrix
## =============================================================================

image2D.matrix <- function (z, x = seq(0, 1, length.out = nrow(z)), 
                   y = seq(0, 1, length.out = ncol(z)), colvar = z, ..., 
                   col = NULL, NAcol = "white", breaks = NULL,
                   border = NA, facets = TRUE, contour = FALSE, 
                   colkey = NULL, resfac = 1, clab = NULL, 
                   lighting = FALSE, shade = NA, ltheta = -135, lphi = 0,
                   theta = 0, rasterImage = FALSE,
                   add = FALSE, plot = TRUE) {

  if (rasterImage & theta != 0)
    stop ("cannot combine 'rasterImage' and 'theta' != 0")
    
  plist <- initplist(add)

  plist <- add2Dplist(plist, "image", z = z, x = x, y = y, 
            col = col, NAcol = NAcol, breaks = breaks, border = border,
            facets = facets, contour = contour, colkey = colkey, resfac = resfac,
            clab = clab, theta = theta, rasterImage = rasterImage,
                    ...)
  setplist(plist)
  if (!plot) return()
  
  if (is.character(z)) {
    ImageNULL (z = NULL, x = x, y = y, ..., col = z, NAcol = NAcol, 
              border = border, facets = facets,
              rasterImage = rasterImage, angle = theta, add = add)   
   return(invisible())
  } 

  if (is.null(col) & is.null(breaks))
    col <- jet.col(100)
  else if (is.null(col))
    col <- jet.col(length(breaks)-1)

  if (length(col) == 1)
    if (is.na(col)) 
      col <- NULL
  
  dots <- splitpardots(list(...))
  dotimage <- dots$main
  dotother <- dots$points

  iscolkey <- is.colkey(colkey, col)       
  if (iscolkey) {
    colkey <- check.colkey(colkey)
    if (! add)       
      plist$plt$main <- colkey$parplt
    setplist(plist)
    colkey$breaks <- breaks
  }
  par (plt = plist$plt$main)


  breaks <- check.breaks(breaks, col)

 # check contours
  iscontour <- ! is.null(contour)
  if (length(contour) == 0) 
    iscontour <- FALSE
  else if (is.null(names(contour)) & is.logical(contour[[1]][1]))
    if (contour[[1]][1] == FALSE) 
      iscontour <- FALSE
  else if (! is.list(contour)) 
    contour <- list()   

 # x- and y-values
  if (length(dim(x)) > 2 | length(dim(y)) > 2)
    stop("'x' or 'y' cannot be an array")
    
  if (is.null (x)) 
    x <- seq(0, 1, length.out = nrow(z))

  if (is.null (y)) 
    y <- seq(0, 1, length.out = ncol(z))

  if (is.matrix(x)) 
    if (! is.matrix(y)) 
      y <- matrix(nrow = nrow(x), ncol = ncol(x), data = y, byrow = TRUE)
     
  if (is.matrix(y)) 
    if (! is.matrix(x)) 
      x <- matrix(nrow = nrow(y), ncol = ncol(y), data = x)

  if (!lighting & is.na(shade)) 
    colvar <- NULL
  else if (any(dim(colvar) - dim(z)) != 0)
    stop("'colvar' and 'z' not compatible")  
 # change resolution
  if (any(resfac != 1)) { 
   if (lighting | !is.na(shade)) 
     res <- changeres(resfac, x, y, z, colvar)
   else
     res <- changeres(resfac, x, y, z)
    x <- res$x
    y <- res$y
    z <- res$z
   if (lighting | !is.na(shade)) 
     colvar <- res$colvar
  }
 
  if (iscontour) {
    if (is.matrix(x))
      stop ("cannot add contour if 'x' or 'y' is a matrix")
    contour$x <- x
    contour$y <- y
    if (dots$clog)
      contour$drawlabels = FALSE        # to avoid strange values
  }
 # rotate 
  rotate <- FALSE
  if (theta != 0 & ! rasterImage) {        
    if (is.vector(x)) {
      rotate <- TRUE
      x <- matrix (nrow = nrow(z), ncol = ncol(z), data = x)
      y <- matrix (nrow = nrow(z), ncol = ncol(z), data = y, byrow = TRUE)
    }

    th <- theta*pi/180
    Mat <- matrix(nrow = 2, data = c(cos(th), sin(th), -sin(th), cos(th)))
    XY <- cbind(as.vector(x), as.vector(y))%*%Mat
    x <- matrix(nrow = nrow(z), ncol = ncol(z), data = XY[, 1])
    y <- matrix(nrow = nrow(z), ncol = ncol(z), data = XY[, 2])
  }

 # Check for decreasing values of x and y    
  if (! is.matrix(x) & all(diff(x) < 0)) {     
    if (is.null(dotimage$xlim)) 
      dotimage$xlim <- rev(range(x))
    x <- rev(x)
    z <- z[nrow(z):1, ]
    if (! is.null(colvar))
      colvar <- colvar[nrow(colvar):1, ]
    
    if (iscontour)
      contour$x <- x
  }
  
  if (! is.matrix(y) & all(diff(y) < 0)) {    
    if (is.null(dotimage$ylim)) 
      dotimage$ylim <- rev(range(y))
    y <- rev(y)
    z <- z[, (ncol(z):1)]
    if (! is.null(colvar))
      colvar <- colvar[, (ncol(colvar):1)]
    if (iscontour)
      contour$y <- y
   }

  useimage <- TRUE     # default it to use the image function
  lightshade <- FALSE
  height <- NULL
  if (lighting | !is.na(shade)) {
    if (! rasterImage & ! is.matrix(x)) {
      xy <- mesh(x,y)
      x <- xy$x
      y <- xy$y  
    }
    if (! rasterImage)
      useimage <- TRUE

    lightshade <- TRUE
    height <- z
    z <- colvar
  }
  Extend <- TRUE

  if (is.matrix(x)) {
    if (any (dim(x) - dim(y) != 0))
      stop("matrices 'x' and 'y' not of same dimension") 
    if (any (dim(x) - dim(z) > 0)) {
      if ((nrow(x) - nrow(z)) != (ncol(x) - ncol(z)))
        stop("matrices 'x' or 'y' and 'z' not compatible - should either be of dim(z) or dim(z)+1")
      if (any(dim(x) - dim(z) > 1))
        stop("matrices 'x' or 'y' and 'z' not compatible - should either be of dim(z) or dim(z)+1")
      Extend <- FALSE
    }
    useimage <- FALSE
  } 
  
  
 # log transformation of z-values (can be specified with log = "c", or log = "z"
  zlog <- FALSE 
  if (! is.null(dots$clog)) 
    zlog <- dots$clog  
  if (! is.null(dotimage$log)) {
    if (length(grep("z", dotimage[["log"]])) > 0) {
      dotimage[["log"]] <- gsub("z", "", dotimage[["log"]])
      zlog <- TRUE
    }
    if (length(grep("c", dotimage[["log"]])) > 0) {
      dotimage[["log"]] <- gsub("c", "", dotimage[["log"]])
      zlog <- TRUE
    }
    if (dotimage[["log"]] == "")
      dotimage[["log"]] <- NULL
  }
  if (zlog) 
    z <- log(z)

 # labels
  if (is.null(dotimage[["xlab"]])) 
    dotimage[["xlab"]] <- "x"
  if (is.null(dotimage[["ylab"]])) 
    dotimage[["ylab"]] <- "y"

 # z ranges
  zlim <- dotimage[["zlim"]]
  if (is.null(zlim)) 
    zlim <- dotother[["clim"]]
  else if (!is.null(dotother[["clim"]])) 
    stop ("only one of 'zlim' and 'clim' can be specified")
    
  dotimage[["zlim"]] <- dotother[["clim"]] <- NULL
  
  if (is.null(zlim)) {
  
    if (length(which(!is.na(z))) == 0)
      zlim <- c(0, 1)
    else if (is.null(breaks))
      zlim <- range(z, na.rm = TRUE)
    else
      zlim <- range(breaks, na.rm = TRUE)
  } else
    if (zlog) 
      zlim  <- log(zlim )

  if (! is.null(dots$alpha))
    col <- setalpha(col, dots$alpha)
  colkeyZlim <- zlim
  colkeyCol  <- col
  


 # Colors for values = NA 
  if (! is.null(NAcol) ) {             #    any (is.na(z)) &
    if (! is.null(breaks)) {
      col <- c(NAcol, col, NAcol)
      breaks <- c(min(c(z, breaks), na.rm = TRUE)-1, breaks,
                  max(c(z, breaks), na.rm = TRUE)+1)
       z[z < min(zlim)] <- NA
       z[z > max(zlim)] <- NA
       z[is.na(z)] <- breaks[1]

    } else {
      nc <- length(col)
      CC <- checkcolors(z, col, NAcol, zlim)
      zlim  <- CC$lim
      col <- CC$col
      z <- CC$colvar
    }
  }

  if (! facets | is.na(facets)) {
    useimage <- FALSE
    if (! is.matrix(x)) {
      xy <- mesh(x,y)
      x <- xy$x
      y <- xy$y  
    }
  }

  if (! useimage | rasterImage) {  # use colored polygons if x- and y are matrix

    # create colors
    Col    <- variablecol (z, col, NAcol, zlim, breaks)

    if (! is.na(shade) | lighting) 
      Col <- facetcolsImage(x, y, height, dotimage[["xlim"]], dotimage[["ylim"]], 
        NULL, shade, lighting, dots$alpha, ltheta, lphi, Col, NAcol)
        
 # empty plot
    dotimage$type <- "n"
    dotimage$xaxs <- "i"
    dotimage$yaxs <- "i"
    if (!add) {
      do.call("plot", c(alist(x = range(x), y = range(y)), dotimage))
    
# This used to make the background = NAcol - removed...
#      plotrect <- !is.null(NAcol)
#      if (plotrect) 
#        if (NAcol != "white") {
#          usr <- par("usr") 
#          rect(usr[1], usr[3], usr[2], usr[4], col = NAcol)
#        }    
    }
  }
  
  if (! useimage ) {
    # function to draw polygon
    poly <- polyfill2D (x, y, Col, facets, border, dots$lwd, dots$lty, Extend)

    dots$lwd <- NULL
    do.call("polygon", c(alist(poly$x, poly$y, lwd = poly$lwd, 
      border = poly$border, col = poly$col), 
                       dotother))
 
  } else if (rasterImage) {
    Col <- matrix(nrow = nrow(z), data = Col)
    addraster (x, y, Col, dotimage[["xlim"]], dotimage[["ylim"]],
      theta, dotother)

  } else {
    dotimage$breaks <- breaks
    do.call("image", c(alist(z = z, x = x, y = y, col = col, add = add,
        zlim = zlim), dotimage))
  }

  if (useimage & !is.na(border)) {
    do.call("abline", c(alist(h = 0.5*(y[-1]+y[-length(y)]), col = border), dotother))
    do.call("abline", c(alist(v = 0.5*(x[-1]+x[-length(x)]), col = border), dotother))
  }
  if (is.null(dotimage$frame.plot)) {
    if (!add)
      box()
  } else if (dotimage$frame.plot)
    box()
    
  # contours
  if (iscontour) {
    if (zlog) 
      if (!is.null(contour$levels)) 
        contour$levels <- log(contour$levels)
    if (! is.null(contour$col) &! is.null(contour$alpha)) 
      contour$col <- setalpha(contour$col, contour$alpha)
  
    if (! rotate)
      do.call("contour", c(list(z = z, add = TRUE), contour))
    else {    # first calculate contours on unrotated values, then transform
      line.list <- do.call("contourLines", c(alist(z), contour))
      templines <- function(clines) 
         lines(cbind(clines[[2]], clines[[3]])%*%Mat)              
        invisible(lapply(line.list, templines))
    }
  }  
  
  if (iscolkey)  {
    drawcolkey(colkey, colkeyCol, colkeyZlim, clab, zlog) 
    par(plt = plist$plt$ori)  
  }               
  par(mar = par("mar")) 
   
}
## =============================================================================
## add rasterImage to plot
## =============================================================================

addraster <- function (x, y, col, xlim, ylim, angle, dots) { 
  
  if (is.matrix(x) | is.matrix(y))
    stop("'x' or 'y' cannot be a matrix if rasterImage is used")

# check the x- and y-values, to be ~equally spaced and monotonously increasing/decreasing
  dx <- diff(x)
  if (any(dx == 0) | max(sign(dx)) != min(sign(dx)))
    stop("'x'-values should be increasing or decreasing, not constant")
    
  if (max(dx)/min(dx) > 1.1)
    stop("'x' should be quasi-equally spaced if  rasterImage is used")
  dy <- diff(y)
  if (any(dy == 0) | max(sign(dy)) != min(sign(dy)))
    stop("'y'-values should be increasing or decreasing, not constant")
  if (max(dy)/min(dy) > 1.1)
    stop("'y' should be quasi-equally spaced if  rasterImage is used")

  if (is.null(xlim)) xlim <- range(x)
  if (is.null(ylim)) ylim <- range(y)

  if (sign(dx[1]) != sign(diff(xlim)))  
    col <- col[nrow(col):1, ]
  if (sign(dy[1]) == sign(diff(ylim)))  
    col <- col[, ncol(col):1]
  col <- t(col)

  if (sign(diff(xlim)) == 1)  
    xlim <- range(x)
  else
    xlim <- rev(range(x))
     
  if (sign(diff(ylim)) == 1)  
    ylim <- range(y)
  else
    ylim <- rev(range(y))
    
  do.call("rasterImage", c(alist(as.raster(col), xlim[1], 
    ylim[1], xlim[2], ylim[2], angle = angle), dots))
}

## =============================================================================
## image2D function, z = NULL, col is a matrix of colors
## =============================================================================

ImageNULL <- function(z = NULL,
                      x = seq(0, 1, length.out = nrow(col)),
                      y = seq(0, 1, length.out = ncol(col)), ...,
                      col, NAcol = "white",
                      border = NA, facets = TRUE,
                      rasterImage = FALSE, angle, add) {

  # check colors
  if (! is.character(col) | ! is.matrix(col))
    stop ("'col' should be a matrix of colors if 'z' is NULL")

 # The plotting arguments
  dots <- splitpardots(list(...))
  dotimage <- dots$main
  dotother <- dots$points

  if (! is.null(dots$alpha)) {
    DD <- dim(col)
    col <- matrix (nrow = DD[1], data = alpha.col(col, dots$alpha))
  }

 # x- and y-values
  if (length(dim(x)) > 2 | length(dim(y)) > 2)
    stop("'x' or 'y' cannot be an array")

  Nr <- nrow(col)
  Nc <- ncol(col)
  if (is.null (x))
    x <- seq(0, 1, length.out = Nr)

  if (is.null (y))
    y <- seq(0, 1, length.out = Nc)

 # labels
  if (is.null(dotimage[["xlab"]]))
    dotimage[["xlab"]] <- "x"
  if (is.null(dotimage[["ylab"]]))
    dotimage[["ylab"]] <- "y"

 # Colors for values = NA
  col[is.na(col)] <- NAcol

 # empty plot
  dotimage$type <- "n"
  dotimage$xaxs <- "i"
  dotimage$yaxs <- "i"
  if (!add) {
    do.call("plot", c(alist(x = range(x), y = range(y)), dotimage))

    plotrect <- !is.null(NAcol)
    if (plotrect)
      if (NAcol != "white") {
        usr <- par("usr")
        rect(usr[1], usr[3], usr[2], usr[4], col = NAcol)
      }
  }
                                       
  if (! rasterImage) {
    if (! is.matrix(x))
      x <- matrix(nrow = Nr, ncol = Nc, data = x)
    if (! is.matrix(y))
      y <- matrix(nrow = Nr, ncol = Nc, data = y, byrow = TRUE)

    if (any (dim(x) - dim(y) != 0))
      stop("matrices 'x' and 'y' not of same dimension")
    if (any (dim(x) - dim(col) != 0))
      stop("matrices 'x' or 'y' and 'col' not of same dimension")

  # function to draw polygon
    poly <- polyfill2D(x, y, col, facets, border, dots$lwd, dots$lty)

    do.call("polygon", c(alist(poly$x, poly$y, lwd = poly$lwd,
      border = poly$border, col = poly$col), dotother))

  } else 
    addraster (x, y, col, dotimage[["xlim"]], dotimage[["ylim"]], 
      angle, dotother)

  if (!add)
    box()

}
## =============================================================================
## image2D function, input is an array
## =============================================================================

image2D.array <- function (z, margin = c(1, 2), subset, ask = NULL, ...) {
  
  DD <- dim(z)
  if (length(DD) != 3)
    stop ("Can only make image of 3-D array, 'z' has dimension ", length(DD))

  if (length(margin) != 2)
    stop ("'margin' should contain two numbers, the x, y subscripts of which to make images")
   
  if ( max(margin) > 3 | min (margin) < 1)
    stop ("indexes in 'margin' should be inbetween 1 and 3")

  index <- (1:3) [- margin]

  if (index > 3 || index <1)
    stop ("'index' to loop over should be inbetween 1 and 3")
  
  x <- 1:DD[index]
  
  if (!missing(subset)){
    if (is.numeric(subset)) { 
      isub <- subset
    } else {
      e <- substitute(subset)
      r <- eval(e, as.data.frame(z), parent.frame())
      if (!is.logical(r))
          stop("'subset' must evaluate to logical")
      isub <- r & !is.na(r)
      isub <- which(isub)
      if (length(isub) == 0)
        stop("cannot continue: nothing selected - check 'subset'")
    }   
  } else isub <- x

  np     <- length(isub)
  ldots  <- list(...)
  
  ## Set par mfrow and ask
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  if (is.null(ldots$main)) {
    title <- names(z)[index][isub]
    if (is.null(title)) 
      title <- isub
  } else 
    title <- rep(ldots$main, length.out = length(isub))

  # outer margin text
  Mtext <- ldots$mtext
  ldots$mtext <- NULL

  i1 <- 1
  for (i in isub) {
    if (index == 1) 
      zz <- z[i, , ] 
    else if (index == 2)
      zz <- z[ ,i , ] 
    else
      zz <- z[ ,, i ]
    if (margin[2] < margin[1])
      zz <- t(zz)

    LL <- c(list(z = zz), ldots)
    LL$main <- title[i1]
    i1 <- i1+1
    do.call(image2D, LL)
  }
  if (! is.null(Mtext))
    mtext(text = Mtext, side = 3, outer = TRUE, line = par("oma")[3]-1 )
  
}

## =============================================================================
## image2D of a list of matrices or arrays
## =============================================================================

image2D.list <- function (z, ...) {
  
# check z: list with similar matrices or arrays of dimension at most 3
  if ( all(c("x", "y", "z") %in% names(z)))  {
    image2D.matrix(z = z$z, x = z$x, y = z$y, ...)
  } else {
    nz     <- length(z)
    classz <- class(z[[1]])

    if (! classz %in% c("matrix", "array"))
      stop ("'z' should be a list with either matrices or arrays")
    
    DD <- dim(z[[1]])
    if (length(DD) > 3 | length(DD) < 2)
      stop ("Can only make image of 2-D or 3-D array, 'z' has dimension ", length(DD))

    for (i in 2 : nz)
      if (any(dim(z[[i]]) - DD != 0))
        stop("elements of 'z' should have the same dimension, check element", i)
  
# Set the mfrow argument - different from the usual
    if ("matrix" %in% classz)  {
      nc <- min(ceiling(sqrt(nz)), 3)
      nr <- min(ceiling(nz/nc), 3)
    } else { # differs from default in that it is not limited to 3
      nc <- ceiling(sqrt(nz))
      nr <- ceiling(nz/nc)
    }
    mfrow <- c(nr, nc)
    par(mfrow = mfrow)

# Plotting arguments
    Ldots <- list(...) 
    Ldots$mfrow <- mfrow
   
    if (!is.null(Ldots$main)) {
      main <- rep(Ldots$main, length.out = nz)
      Ldots$main <- NULL
    } else {
      main <- names(z)
      if (is.null(main)) main <- 1:nz
    }
    ask <- Ldots$ask
    if (is.null(ask)) ask <- TRUE
    Ldots$ask <- NULL

  # ylim and xlim can be lists and are at least two values
    yylim  <- expanddotslist(Ldots$ylim, nz)
    xxlim  <- expanddotslist(Ldots$xlim, nz)
    zzlim  <- expanddotslist(Ldots$zlim, nz)
    zzlab  <- expanddotslist(Ldots$clab, nz)
   
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

 # Display the images
    if ("matrix" %in% classz)  {
     # outer margin text
      Mtext <- Ldots$mtext
      Ldots$mtext <- NULL

      for (i in 1:nz) {
        Ldots$main <- main[i]
        Ldots$xlim <- xxlim[[i]]
        Ldots$ylim <- yylim[[i]]
        Ldots$zlim <- zzlim[[i]]
        Ldots$clab <- zzlab[[i]]
      
        LL <- c(list(z = z[[i]]), Ldots)
        do.call(image2D, LL)
      }
  
      if (! is.null(Mtext))
        mtext(text = Mtext, side = 3, outer = TRUE, line = par("oma")[3]-1 )
   
    } else {  # array
      margin <- Ldots$margin
      Ldots$margin <- NULL
     
      if (is.null(margin)) margin <- 1:2
      if (length(margin) != 2)
        stop ("'margin' should contain two numbers, the x, y subscripts with which to make images")
      if ( max(margin) > 3 | min (margin) < 1)
        stop ("indexes in 'margin' should be inbetween 1 and 3")
      index <- (1:3) [- margin]

      subset <- Ldots$subset
      Ldots$subset <- NULL
      if (!is.null(subset)){
        if (is.numeric(subset)) { 
          isub <- subset
        } else {        e <- substitute(subset)
          r <- eval(e, as.data.frame(z), parent.frame())
          if (!is.logical(r))
            stop("'subset' must evaluate to logical")
          isub <- r & !is.na(r)
          isub <- which(isub)
          if (length(isub) == 0)
            stop("cannot continue: nothing selected - check 'subset'")
        } 
      } else 
        isub <- 1:DD[index]
     
      nisub     <- length(isub)

   # number of empty plots 
      noplot <- prod(mfrow) - nz
      if (noplot == 0) 
        noplot <- NULL 
      else 
        noplot <- 1:noplot

     # outer margin text
      Mtext <- Ldots$mtext
      Ldots$mtext <- NULL
     
      if (! is.null(Mtext)) 
        Mtext <- rep(Mtext, length.out = nisub)
      else
        Mtext <- isub
      pline <- par("oma")[3]-1  
   # loop first over margin, then over data sets
      for (jj in 1:nisub) {
        j <- isub[jj]      
        for (i in 1:nz) {
          if (index == 1) 
            zz <- z[[i]][j, , ] 
          else if (index == 2)
            zz <- z[[i]][ ,j , ] 
          else
            zz <- z[[i]][ ,, j ]
          if (margin[2] < margin[1])
            zz <- t(zz)
       
          Ldots$main <- main[i]
          Ldots$xlim <- xxlim[[i]]
          Ldots$ylim <- yylim[[i]]
          Ldots$zlim <- zzlim[[i]]
          Ldots$clab <- zzlab[[i]]
          LL <- c(list(z = zz), Ldots)
          do.call(image2D, LL)
        }
     # to make sure all figures are drawn
        for (i in noplot) 
          plot(0, type = "n", xlab = "", ylab = "", axes = FALSE, 
              frame.plot = FALSE)
        mtext(text = Mtext[jj], side = 3, outer = TRUE, line = pline)
      }  
    } 
  }
}

## =============================================================================
## Checking and expanding arguments in dots (...) with default
## =============================================================================

expanddots <- function (dots, default, n) {
  dots <- if (is.null(dots)) default else dots
  rep(dots, length.out = n)
}

# lists: e.g. xlim and ylim....
expanddotslist <- function (dots, n) {
  if (is.null(dots)) return(dots)
  dd <- if (!is.list(dots )) list(dots) else dots
  rep(dd, length.out = n)
}


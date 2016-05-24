## =============================================================================
## function that makes a box type in 2D
## =============================================================================

start2Dplot <- function(dots, x, y) {
  dd <- dots
  dd$type <- "n"

  bty <- dots$bty
  dots$bty <- NULL

  if (is.null(bty))
    bty <- "o"

  if (bty %in% c("b2", "g", "bl")) 
    dd$bty <- NULL
    
  if (is.null(dd$xlab))
    dd$xlab <- "x"
      
  if (is.null(dd$ylab))
    dd$ylab <- "y"

  do.call("plot", c(alist(x, y), dd))

  if (bty == "b2")
    grid(col = "grey", lty = 1, lwd = 2)
  else if (bty == "g") {
    pu <- par("usr")
    rect(pu[1], pu[3], pu[2], pu[4],
      col = grey(0.925), border = grey(0.925))
    grid(col = "white", lty = 1, lwd = 2)
  } else if (bty %in% c("bl","bl2")) {
    pu <- par("usr")
    rect(pu[1], pu[3], pu[2], pu[4], col = "black")
    if (bty == "bl2")
      grid(col = "grey", lty = 1, lwd = 2)
  }
  return(dots)
}

## =============================================================================

plot2Dplist <- function (plist, ...) {
  
  checkdots <-  function(pdots, dots, add) {
    if (! add) {
      if (! is.null(dots$xlim)) 
        pdots$xlim <- dots$xlim
      if (! is.null(dots$ylim)) 
        pdots$ylim <- dots$ylim
      if (! is.null(dots$alpha)) 
        pdots$alpha <- dots$alpha
    }
    pdots
  }
  
  img2Dnr <- cont2Dnr <- scat2Dnr <- arr2Dnr <- segm2Dnr <- rect2Dnr <- poly2Dnr <- text2Dnr <- 0
  add <- FALSE
  dots <- list(...)
  p <- plist$twoD
  for (i in 1:length(p$order)) {
    plt <- p$order[i]
    if (plt  == "image") {
      img2Dnr <- img2Dnr + 1
      Dots <- checkdots(p$img2D[[img2Dnr]], dots, add)
      do.call ("image2D", c(alist(add = add), Dots))
    } else if (plt  == "contour") {
      cont2Dnr <- cont2Dnr + 1
      Dots <- checkdots(p$cont2D[[cont2Dnr]], dots, add)
      do.call ("contour2D", c(alist(add = add), Dots))
    } else if (plt == "scatter") {
      scat2Dnr <- scat2Dnr + 1
      Dots <- checkdots(p$scat2D[[scat2Dnr]], dots, add)
      do.call ("scatter2D", c(alist(add = add), Dots))
    } else if (plt %in% c("arrows", "ArrType")) {
      arr2Dnr <- arr2Dnr + 1
      Dots <- checkdots(p$arr2D[[arr2Dnr]], dots, add)
      do.call ("arrows2D", c(alist(add = add), Dots))
    } else if (plt == "segments") {
      segm2Dnr <- segm2Dnr + 1
      Dots <- checkdots(p$segm2D[[segm2Dnr]], dots, add)
      do.call ("segments2D", c(alist(add = add), Dots))
    } else if (plt == "rect") {
      rect2Dnr <- rect2Dnr + 1
      Dots <- checkdots(p$rect2D[[rect2Dnr]], dots, add)
      do.call ("rect2D", c(alist(add = add), Dots))
    } else if (plt == "polygon") {
      poly2Dnr <- poly2Dnr + 1
      Dots <- checkdots(p$poly2D[[poly2Dnr]], dots, add)
      do.call ("polygon2D", c(alist(add = add), Dots))
    } else if (plt  == "text") {
      text2Dnr <- text2Dnr + 1
      Dots <- checkdots(p$text2D[[text2Dnr]], dots, add)
      do.call ("text2D", c(alist(add = add), Dots))
    } 
    add <- TRUE
  }
  invisible(plist)
}

## =============================================================================

add2Dplist <- function(plist, method, ...) {
  dots <- list(...)
  
  if (is.null(plist) | length(plist) == 0) {
    setlim <- c(!is.null(dots$xlim), !is.null(dots$ylim), !is.null(dots$zlim))
    plist <- list(type = "2D", 
      xlim = dots$xlim, ylim = dots$ylim, zlim = dots$zlim, setlim = setlim,
      twoD = list(order = NULL), plt = list(main = par("plt"), ori = par("plt")))
  } 
  
  if (plist$type == "3D")
    plist$type <- "23D" 
#    stop ("cannot merge 2D and 3D plotting functions")
  p <- plist$twoD 
  p$order <- c(p$order, method)
  
  if (method == "image") {
    if (is.null(p$img2Dnr)) {
        p$img2Dnr <- 0
        p$img2D <- list()
    } 
    p$img2Dnr <- p$img2Dnr + 1
    p$img2D[[p$img2Dnr]] <- dots
  }
    
  else if (method == "contour") {
    if (is.null(p$cont2Dnr)) {
        p$cont2Dnr <- 0
        p$cont2D <- list()
    } 
    p$cont2Dnr <- p$cont2Dnr + 1
    p$cont2D[[p$cont2Dnr]] <- dots
  }

  else if (method == "scatter") {
    if (is.null(p$scat2Dnr)) {
        p$scat2Dnr <- 0
        p$scat2D <- list()
    } 
    p$scat2Dnr <- p$scat2Dnr + 1
    p$scat2D[[p$scat2Dnr]] <- dots
  }
    
  else if (method %in% c("arrows", "ArrType")) {
    if (is.null(p$arr2Dnr)) {
        p$arr2Dnr <- 0
        p$arr2D <- list()
    } 
    p$arr2Dnr <- p$arr2Dnr + 1
    p$arr2D[[p$arr2Dnr]] <- dots
  }
    
  else if (method == "segments") {
    if (is.null(p$segm2Dnr)) {
        p$segm2Dnr <- 0
        p$segm2D <- list()
    } 
    p$segm2Dnr <- p$segm2Dnr + 1
    p$segm2D[[p$segm2Dnr]] <- dots
  }

  else if (method == "rect") {
    if (is.null(p$rect2Dnr)) {
        p$rect2Dnr <- 0
        p$rect2D <- list()
    } 
    p$rect2Dnr <- p$rect2Dnr + 1
    p$rect2D[[p$rect2Dnr]] <- dots
  }
  else if (method == "polygon") {
    if (is.null(p$poly2Dnr)) {
        p$poly2Dnr <- 0
        p$poly2D <- list()
    } 
    p$poly2Dnr <- p$poly2Dnr + 1
    p$poly2D[[p$poly2Dnr]] <- dots
  }
  else if (method == "text") {
    if (is.null(p$text2Dnr)) {
        p$text2Dnr <- 0
        p$text2D <- list()
    } 
    p$text2Dnr <- p$text2Dnr + 1
    p$text2D[[p$text2Dnr]] <- dots
  }

  plist$twoD <- p
  plist$colkeyargs <- dots$colkey
  class(plist) <- c("plist","list")
  plist
}

## =============================================================================
## plots (2-D)
## =============================================================================

# x, y, colvar: vector or matrix of same dimension

plot2D <- function(x0, y0, x1, y1, ..., colvar = NULL,
                    col = NULL, NAcol = "white", breaks = NULL,
                    colkey = NULL,
                    clim = NULL, clab = NULL, add = FALSE,
                    plot = TRUE, method = "arrows") {

  plist <- initplist(add)

  plist <- add2Dplist(plist, method, x0 = x0, y0 = y0, x1 = x1, y1 = y1, 
    colvar = colvar, col = col, NAcol = NAcol, breaks = breaks,
    colkey = colkey, clim = clim, clab = clab, ...)
  setplist(plist)

  if (!plot) return()
  dots <- splitpardots(list(...))

 # colors
  breaks <- check.breaks(breaks, col)
  if (! is.null(colvar)) {
    if (is.null(col))
      col <- jet.col(100)

    if (dots$clog) {
      colvar <- log(colvar)
      if (! is.null(clim)) 
        clim <- log(clim)
    }

    iscolkey <- is.colkey(colkey, col)

    if (iscolkey) {
      colkey <- check.colkey(colkey)
      if (! add)       
        plist$plt$main <- colkey$parplt
      setplist(plist)
      colkey$breaks <- breaks
    }

    if (length(colvar) != length(x0))
      stop ("length of 'colvar' should be equal to length of 'x0', 'x1', 'y0' and 'y1'")

    if (is.null(clim))
      clim <- range(colvar, na.rm = TRUE)

    if (! is.null(dots$alpha)) 
      col <- setalpha(col, dots$alpha)
    Col <- variablecol(colvar, col, NAcol, clim, breaks)

  } else  {  # no colvar
    Col <- col
    if (is.null(Col)) Col <- "black"
    if (! is.null(dots$alpha)) 
      Col <- setalpha(Col, dots$alpha)
    iscolkey <- FALSE
  }
  par (plt = plist$plt$main)

  if (! add)
    dots$main <- start2Dplot(dots$main, c(x0, x1), c(y0, y1))
    
  do.call(method, c(alist(x0, y0, x1, y1, col = Col), dots$points))

  if (iscolkey) {
    drawcolkey(colkey, col, clim, clab, dots$clog)
    par(plt = plist$plt$ori)  
  }
  par(mar = par("mar"))

}

## =============================================================================
## specific plot functions (2-D)
## =============================================================================
arrows2D <- function(x0, y0, x1 = x0, y1 = y0,..., colvar = NULL,
                    col = NULL, NAcol = "white", breaks = NULL,
                    colkey = NULL,
                    clim = NULL, clab = NULL, 
                    type = "triangle", add = FALSE, plot = TRUE)  
  plot2D (x0, y0, x1, y1, ..., colvar = colvar,
                    col = col, NAcol = NAcol, breaks = breaks,
                    colkey = colkey, type = type,
                    clim = clim, clab = clab, add = add,
                    plot = plot, method = "ArrType")

## =============================================================================
segments2D <- function(x0, y0, x1 = x0, y1 = y0, ..., colvar = NULL,
                    col = NULL, NAcol = "white", breaks = NULL,
                    colkey = NULL,
                    clim = NULL, clab = NULL, add = FALSE, plot = TRUE) 
  plot2D (x0, y0, x1, y1, ..., colvar = colvar,
                    col = col, NAcol = NAcol, breaks = breaks,
                    colkey = colkey,
                    clim = clim, clab = clab, add = add,
                    plot = plot, method = "segments")

## =============================================================================
rect2D <- function(x0, y0, x1 = x0, y1 = y0, ..., colvar = NULL,
                    col = NULL, NAcol = "white", breaks = NULL,
                    colkey = NULL,
                    clim = NULL, clab = NULL, add = FALSE, plot = TRUE) 
  plot2D (x0, y0, x1, y1, ..., colvar = colvar,
                    col = col, NAcol = NAcol, breaks = breaks,
                    colkey = colkey,
                    clim = clim, clab = clab, add = add,
                    plot = plot, method = "rect")


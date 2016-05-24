
## =============================================================================
## =============================================================================
## Color key functions
## =============================================================================
## =============================================================================

## =============================================================================
## Check if necessary to draw a color key
## =============================================================================

is.colkey <- function(colkey, col) {
  
  if (is.logical(colkey))
    return(colkey)
    
  if (is.list(colkey))
    return(TRUE)
    
  if (! is.null(colkey))
    stop("'colkey' should be a list, a logical or NULL")
 
  iscolkey <- ispresent(col) 
  
  if (iscolkey) {
    if (length(col) == 1) 
      iscolkey <- FALSE
    else if (length(col) == 2 & col[1] == col[2]) 
      iscolkey <- FALSE
  }
  
  return(iscolkey)
}

## =============================================================================
## function to extract default parameter values if not overruled
## =============================================================================

overrulepar <- function(main, subset) {
  nmsC <- names(main)
  main[(namc <- names(subset))] <- subset
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
     warning("unknown names in colkey parameter subset: ", paste(noNms, 
            collapse = ", "))
  return(main)
}

## =============================================================================
## color key parameter check
## =============================================================================

check.colkey <- function(colkeypar, add = FALSE) {    
   
  if (!is.list(colkeypar))
    colkeypar <- list()
      
  parameter <- list(side = 4, plot = TRUE,
      length = 1, width = 1, dist = 0, shift = 0, addlines = FALSE,
      col.clab = NULL, cex.clab = par("cex.lab"), 
      side.clab = NULL, line.clab = NULL, adj.clab = NULL, font.clab = NULL, 
      at = NULL, labels = TRUE, tick = TRUE, line = NA, 
      pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, 
      lwd.ticks = 1, col.box = NULL, col.axis = NULL, 
      col.ticks = NULL, hadj = NA, padj = NA, 
      cex.axis = par("cex.axis"), mgp = NULL, 
      tck = NULL, tcl = NULL, las = NULL)
                             
  colkeypar$parleg <- colkeypar$parplt <- NULL
  
  colkey <- overrulepar(parameter, colkeypar)
  
  if (is.numeric(colkey$labels))
    colkey$labels <- as.logical(colkey$labels)
  if (is.null(colkey$side)) 
    colkey$side <- 4

 # plt parameters of legend
  colkey <- key.parleg(colkey, add)
  return(colkey)
}

## =============================================================================
   
key.parleg <- function(colkey, add) {   # the plotting parameters                                                       

  dw     <- 0.03*colkey$width
  rp     <- colkey$shift
   
  parplt <- par("plt")
    
  if (colkey$side == 1) {
    if (! add) 
      parplt <- par("plt") + c(0, 0, 0.145, 0)
    dd     <- 0.145 + colkey$dist
    dp     <- (parplt[2] - parplt[1]) * (1-colkey$length)/2 
    parleg <- c( parplt[1] + dp+rp, parplt[2] - dp+rp, 
                 parplt[3] - dd, parplt[3] - dd + dw)
  
  } else if (colkey$side == 2) {
    if (! add) 
      parplt <- par("plt") + c(0.1, 0, 0, 0)
    dd <- 0.1  + colkey$dist
    dp <- (parplt[4] - parplt[3]) * (1-colkey$length)/2 
    parleg <- c( parplt[1] - dd, parplt[1] - dd + dw, 
                 parplt[3] + dp+rp, parplt[4] - dp+rp)
  
  } else if (colkey$side == 3) {
    if (! add) 
      parplt <- par("plt") - c(0, 0, 0, 0.08)
    dd <- 0.02 + colkey$dist
    dp <- (parplt[2] - parplt[1]) * (1-colkey$length)/2 
    parleg <- c( parplt[1] + dp+rp,  parplt[2] - dp+rp, 
                 parplt[4] + dd, parplt[4] + dd + dw)
  
  } else if (colkey$side == 4) {
    if (! add) 
      parplt <- par("plt") - c(0, 0.08, 0, 0)
    dd <- 0.02 + colkey$dist
    dp <- (parplt[4] - parplt[3]) * (1-colkey$length)/2
    parleg <- c(parplt[2] + dd, parplt[2] + dd + dw, 
                parplt[3] + dp+rp, parplt[4] - dp+rp)
  }
   
  colkey$parleg <- check.plt(parleg)
  colkey$parplt <- check.plt(parplt)
  return(colkey)
}

## =============================================================================
## function to save the color key settings in the plotting list
## =============================================================================

plistcolkey <- function (plist, colkeypar, col, zlim, zlab = NULL, 
                         zlog = FALSE, New = TRUE,
                         type = "scatter3D", breaks)  {
  if (is.null(plist$colkey))  {
    plist$colkey <- list()
    plist$numkeys <- 1
  } else
    plist$numkeys <- plist$numkeys + 1
  if (! is.null(breaks))
    colkeypar$breaks <- breaks

  plist$colkey[[plist$numkeys]] <- list(par = colkeypar, col = col, 
    clim = zlim, clab = zlab, clog = zlog, New = New, type = type)
  plist
} 
                              
## =============================================================================
## functions to draw the color key
## =============================================================================

drawallcols <- function(plist) {
  for (colkey in plist$colkey) 
    drawcolkey(colkey$par, colkey$col, colkey$clim, colkey$clab, 
                        colkey$clog, colkey$New)
}     

## =============================================================================

drawcolkey <- function (colkeypar, col, clim, clab = NULL, 
                        clog = FALSE, New = TRUE) {     

  if (!colkeypar$plot) return()
  parleg <- check.plt(colkeypar$parleg)

  Plt <- par(plt = parleg)   
  PP  <- par()
  if (New) 
    par(new = TRUE)

  usr <- par("usr")
  col.clab <- colkeypar$col.clab
  cex.clab <- colkeypar$cex.clab
  side.clab <- colkeypar$side.clab
  line.clab <- colkeypar$line.clab
  adj.clab  <- colkeypar$adj.clab
  font.clab <- colkeypar$font.clab
  addlines  <- colkeypar$addlines

  if (is.null(cex.clab))
    cex.clab <- par("cex.lab")
  ix <- 1
  minz <- min(clim)
  maxz <- max(clim)
# the parameters for the axis
  axispar <- colkeypar
  nbins <- length(col)
  binwidth <- (maxz - minz)/nbins
  if (! is.null(axispar$breaks)) {
    nbreaks <- length(axispar$breaks)
    axispar$labels <- axispar$breaks
    minz <- 0.5
    maxz <- nbreaks+0.5
    clim <- c(minz, maxz)
    binwidth <- (maxz - minz)/nbins
    iyb <- seq(minz, maxz, by = binwidth)
    if (clog) iyb <- exp(iyb)

  }
  iy <- IY <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  if (clim[1] > clim[2])
    col <- rev(col)
  if (clog) {
    iy <- exp(iy)
    if (colkeypar$side %in% c(2, 4))  {
      Log <- "y"
    } else {
      Log <- "x"
    }  
  } else Log <- ""
  iz <- matrix(IY, nrow = 1, ncol = length(iy))

  if (! is.numeric(cex.clab)) 
    cex.clab <- 1.

  if (! is.null(axispar$breaks))
    axispar$at <- iyb

# remove arguments not in axis function
  axispar$side <- axispar$length <- axispar$width <- axispar$plot <- NULL
  axispar$parleg <- axispar$parplt <- axispar$dist <- NULL
  axispar$shift <-axispar$col.box <- axispar$breaks <- NULL
  axispar$col.clab <- axispar$cex.clab <- axispar$side.clab <- NULL
  axispar$line.clab <- axispar$adj.clab <- axispar$font.clab <- NULL
  axispar$addlines <- NULL
  
  if (colkeypar$side %in% c(2, 4)) {
    ylim <- clim
    if (Log == "y") ylim <- exp(ylim)

    image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", log = Log,
        ylab = "", col = col, main = "", ylim = ylim)
    if (addlines)  
      abline(h = seq(mean(iy[1:2]), mean(iy[(nbins-1):nbins]), length.out = nbins-1))
    do.call("axis", c(list(side = colkeypar$side, mgp = c(3, 1, 0), las = 2), 
            axispar))
                                                  
  } else {

    xlim <- clim
    if (Log == "x") xlim <- exp(xlim)
    image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", log = Log,
        ylab = "", col = col, main = "", xlim = xlim)
    if (addlines) 
      abline(v = seq(mean(iy[1:2]), mean(iy[(nbins-1):nbins]), length.out = nbins-1))
    do.call("axis", c(list(side = colkeypar$side, mgp = c(3, 1, 0), las = 1), 
         axispar))
  }
  if (is.null(side.clab))
    title(main = clab, cex.main = cex.clab, col.main = col.clab, 
      line = line.clab, adj = adj.clab, font = font.clab)
  else if (side.clab == 1)
    title(xlab = clab, cex.lab = cex.clab, col.lab = col.clab, 
      line = line.clab, adj = adj.clab, font.lab = font.clab)
  else if (side.clab == 2)
    title(ylab = clab, cex.lab = cex.clab, col.lab = col.clab, 
      line = line.clab, adj = adj.clab, font.lab = font.clab)
  else if (side.clab == 3)
    title(main = clab, cex.main = cex.clab, col.main = col.clab, 
      line = line.clab, adj = adj.clab, font = font.clab)
  else {
    if (is.null(adj.clab))
      adj.clab <- NA
    if (is.null(line.clab))
      line.clab <- 2
    mtext(side = side.clab, text = clab, cex = cex.clab, 
        col = col.clab, line = line.clab, adj = adj.clab, font = font.clab)
  }    
  if (clog) {
    if (colkeypar$side %in% c(2, 4))  
      par (ylog = FALSE)
    else 
      par (xlog = FALSE)      
  }

  box(col = colkeypar$col.box) 
  par(plt = Plt)
  par(usr = usr)
  par(xlog = PP$xlog)
  par(ylog = PP$ylog)
  
  if (New) 
    par(new = FALSE)  
}

## =============================================================================
## =============================================================================
## R-function to draw a color key
## =============================================================================
## =============================================================================

colkey <- function(col = NULL, clim, clab = NULL, clog = FALSE, add = FALSE,
                   cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
                   line.clab = NULL, adj.clab = NULL, font.clab = NULL,
                   side = 4, length = 1, width = 1, 
                   dist = 0, shift = 0, addlines = FALSE, 
                   breaks = NULL, at = NULL, labels = TRUE,
                   tick = TRUE, line = NA, pos = NA, outer = FALSE,  
                   font = NA, lty = 1, lwd = 1, lwd.ticks = 1, 
                   col.axis = NULL, col.ticks = NULL, col.box = NULL,
                   hadj = NA, padj = NA, cex.axis = par("cex.axis"),
                   mgp = NULL, tck = NULL, tcl = NULL, las = NULL) {

  if (is.null(col))
    if (is.null(breaks))
      col <- jet.col(100)
    else
      col <- jet.col(length(breaks) - 1)
  breaks <- check.breaks(breaks, col)
  if (! is.null(breaks))
    clim <- range(breaks)
  colkey <- list(side = side, plot = TRUE, length = length, width = width, 
      dist = dist, shift = shift, addlines = addlines, 
      cex.clab = cex.clab, col.clab = col.clab,
      side.clab = side.clab, line.clab = line.clab, adj.clab = adj.clab,
      font.clab = font.clab, breaks = breaks,
      at = at, labels = labels, tick = tick, line = line, 
      pos = pos, outer = outer, font = font, lty = lty, lwd = lwd, 
      lwd.ticks = lwd.ticks, col.box = col.box, col.axis = col.axis, 
      col.ticks = col.ticks, hadj = hadj, 
      padj = padj, cex.axis = cex.axis,
      mgp = mgp, tck = tck, tcl = tcl, las =las)

  if (is.numeric(colkey$labels))
    colkey$labels <- as.logical(colkey$labels)
  if (is.null(colkey$side)) colkey$side <- 4
                                                          
  dw     <- 0.03*colkey$width

  parplt <- par("plt") 
  
  if (! add & colkey$side %in% c(1, 3)) {
      py <- 0.5*(parplt[3] + parplt[4])
      dp <- (parplt[2] - parplt[1]) * (1-colkey$length)/2
      colkey$parleg <- c( parplt[1]+dp, parplt[2]-dp, py - dw/2, py+dw/2)

  } else if (! add & colkey$side %in% c(2, 4)) {
      px <- 0.5*(parplt[1] + parplt[2])
      dp <- (parplt[4] - parplt[3])*(1-colkey$length)/2
      colkey$parleg <- c( px - dw/2, px +dw/2, parplt[3]+dp, parplt[4]-dp)

  } else colkey <- key.parleg(colkey, add = TRUE)
  
  colkey$parplt <- parplt
  
  if (clog) 
    clim <- log(clim) 
  
  drawcolkey (colkey, col, clim = clim, clab = clab, 
                        clog = clog, New = add)

  par(mar = par("mar")) # to prevent R from setting defaultplot = false

}

## =============================================================================
## checks the validity of the plotting arguments "plt"
## =============================================================================

check.plt <- function(plt) {
  if (!(plt[1] < plt[2] & plt[3] < plt[4]))
    stop("figure margins too large")
  eps <- 1e-10
  if (!(plt[1] > -eps & plt[2] < 1+eps & plt[3] > -eps & plt[4] < 1+eps))
    stop("plot region too large")
  return(plt)  
}

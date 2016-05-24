
## =============================================================================
## =============================================================================
## Color key functions
## =============================================================================
## =============================================================================

## =============================================================================
## Is it necessary to draw a color key?
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
 
  iscol <- ispresent(col) 
  
  if (iscol) {
    if (length(col) == 1) 
      iscol <- FALSE
    else if (length(col) == 2 & col[1] == col[2]) 
      iscol <- FALSE
  }
  
  return(iscol)
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
      
  parameter <- list(side = 4, #plot = TRUE,
      length = 1, width = 1, dist = 0, shift = 0, addlines = FALSE,
      col.clab = NULL, cex.clab = par("cex.lab"), 
      side.clab = NULL, line.clab = NULL, adj.clab = NULL, font.clab = NULL, 
      at = NULL, labels = TRUE, tick = TRUE, line = NA, 
      pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, 
      lwd.ticks = 1, col.box = NULL, col.axis = NULL, 
      col.ticks = NULL, hadj = NA, padj = NA, 
      cex.axis = par("cex.axis"), mgp = NULL, 
      tck = NULL, tcl = NULL, las = NULL)
                                   
  colkey <- overrulepar(parameter, colkeypar)
  if (is.numeric(colkey$labels))
    colkey$labels <- as.logical(colkey$labels)
  if (is.null(colkey$side)) 
    colkey$side <- 4

 # impose the distance of the legend
  if (is.null(colkeypar$dist)) {
    colkey$dist <- -0.075     
    if (colkey$side == 3)
      colkey$dist <- 0     
  }
  
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
      parplt <- par("plt") + c(0, 0, 0.15, 0)
    dd     <- 0.15 + colkey$dist
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

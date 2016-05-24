## =============================================================================
## =============================================================================
## General utility functions
## =============================================================================
## =============================================================================

## =============================================================================
## something can be toggled off by putting it = FALSE, NA, NULL
## =============================================================================

ispresent <- function(var) {
  if (length(var) > 1) 
    return(TRUE)

  ispresent <- TRUE
  if (is.null(var))
    ispresent <- FALSE
  else if (is.na(var))
    ispresent <- FALSE
  else if (length(var) == 1)
    if (is.logical(var)) 
      if (!var)
        ispresent <- FALSE
  return(ispresent)
}

## =============================================================================
## List to draw legend/contour?
## =============================================================================

check.args <- function(ll, dots) {
  
  addit <- ! is.null(ll)
 
  if (length(ll) == 0) 
    addit <- FALSE
  
  else if (is.logical(ll[[1]]))
    if (length(ll[[1]]) == 1)      
      if (ll[[1]] == FALSE) 
        addit <- FALSE
  
  side <- 1
  colkey <- FALSE  
  
  if (addit) {     # should have at least side argument
    if (is.list(ll)) {
      if (!is.null(ll$side))
        side <- ll$side
      ll$side <- NULL  
    } else ll <- list() 
    ll <- c(ll, dots)
    if (!is.null(ll$colkey)) {
      colkey <- ll$colkey
      ll$colkey <- NULL
    }  
  }
  list(add = addit, side = side, args = ll, colkey = colkey)
}

## =============================================================================
## Functions that account for occurrence of decreasing values...
## =============================================================================

FindInterval <- function(x, vec, ...) {

  if (all(diff(vec) < 0)) {# swap
    vec <- rev(vec)
    res <- c(length(vec):1) [findInterval(x, vec, ...)]-1
  } else 
    res <- findInterval(x, vec, ...)
  res [ res == 0] <- 1
  res
}

## =============================================================================
## .. and of NAs
## =============================================================================

Approx <- function(x, y, xout, ...) {

  if (all(is.na(x)) | all(is.na(y))) 
    return(list(y = rep(NA, length = length(xout)), x = xout))

  if (diff(range(x, na.rm = TRUE)) == 0) {
    warning("Warning in approx: all 'x' values are the same - returning 'NA'")
    return(list(y = rep(NA, length = length(xout)), x = xout))
  }
  if (any(is.na(c(x, y)))) {
    ii <- unique(c(which(is.na(x)), which(is.na(y))))
    x <- x[-ii]
    y <- y[-ii]
  } 
  approx(x, y, xout, ...)
}

## =============================================================================
## =============================================================================
## Colors
## =============================================================================
## =============================================================================

check.breaks <- function (breaks, col) {
  if (! is.null(breaks)) {
     nbreaks <- length(breaks)
     if (length(col) != nbreaks-1)
       stop("must have one more break than col - suggest to use jet.col(", nbreaks-1, ")")

  if (any(!is.finite(breaks)))
            stop("'breaks' must all be finite")
  if (is.unsorted(breaks)) {
            warning("unsorted 'breaks' will be sorted before use")
            breaks <- sort(breaks)
         }
  }
  return(breaks)
}

## =============================================================================
## Generates color vector based on variable values
## =============================================================================

variablecol <- function(colvar, col, NAcol, clim, breaks) {
 
  if (is.null(breaks)) {
   ncol <- length(col)
   rn <- clim[2] - clim[1]
   ifelse (rn != 0, Col <- col[1 + trunc((colvar - clim[1])/rn *
      (ncol - 1))], Col <- rep(col[1], ncol))
  } else {
      zi <- .bincode(colvar, breaks, TRUE, TRUE)
      Col <- col[zi]
  }

  Col[is.na(Col)] <- NAcol
  return(Col)
}

## =============================================================================
## =============================================================================
## DOTS splitted, expanded, extracted. etc..
## =============================================================================
## =============================================================================

## =============================================================================
## Expanding arguments in dots  (...)
## =============================================================================

repdots <- function(dots, n) {
  if (is.function(dots)) { 
    return(dots)
  } else 
    return(rep(dots, length.out = n))
}

setdots <- function(dots, n) 
  lapply(dots, repdots, n)

## =============================================================================
## Extracting element 'index' from dots  (...)
## =============================================================================

extractdots <- function(dots, index) {
  
  ret <- lapply(dots, "[", index)
  ret <- lapply(ret, unlist) # flatten list
  return(ret)
}

## =============================================================================
## set mfrow and ask
## =============================================================================

setplotpar <- function(ldots, nv, ask) {

  nmdots <- names(ldots) 

  # nv = number of variables to plot
  if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
    nc <- min(ceiling(sqrt(nv)), 3)
    nr <- min(ceiling(nv/nc), 3)
    mfrow <- c(nr, nc)
  } else if ("mfcol" %in% nmdots)
    mfrow <- rev(ldots$mfcol)
  else 
    mfrow <- ldots$mfrow

  if (! is.null(mfrow))  
    mf <- par(mfrow = mfrow)

  ## interactively wait if there are remaining figures
  if (is.null(ask))
    ask <- prod(par("mfrow")) < nv && dev.interactive()

  return(ask)
}

## =============================================================================
## Split plotting parameters in general (main) and point parameters
## =============================================================================

splitpardots <- function(dots) {

  clog <- dots$clog
  if (is.null(clog)) { 
    clog <- FALSE
    if (! is.null(dots$log)) {
      if (length(grep("c", dots[["log"]])) > 0) {
        dots[["log"]] <- gsub("c", "", dots[["log"]])
        if (dots[["log"]] == "")
          dots[["log"]] <- NULL
        clog <- TRUE
      } 
    } 
  }

  nmdots <- names(dots)

  # plotting parameters : split in plot parameters and point parameters
  plotnames <- c("xlab", "ylab", "zlab", "xlim", "ylim", "zlim", 
                 "main", "sub", "log", "asp", "bty", 
                 "ann", "axes", "frame.plot", "panel.first", "panel.last",
                 "cex.lab", "col.lab", "font.lab",
                 "cex.axis", "col.axis", "font.axis", 
                 "cex.main", "col.main", "font.main")

  # plot.default parameters
  ii <- names(dots) %in% plotnames
  dotmain <- dots[ii]

  # point parameters
  ip <- !names(dots) %in% c(plotnames, "add", "clog")
  dotpoints <- dots[ip]
  list (points = dotpoints, main = dotmain, add = dots$add, clog = clog)

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


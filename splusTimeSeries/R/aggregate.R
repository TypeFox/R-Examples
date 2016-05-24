"aggregateSeries" <- 
function(x, pos, FUN, moving = FALSE, together = FALSE,
         drop.empty = TRUE, include.ends = 
	FALSE, adj, offset, colnames, by, k.by = 1, week.align = NULL,
         holidays = 
	timeDate(), align.by = TRUE, incr = 1, ...)
{
  ## Aggregate series object to new positions
  origpos <- positions(x)
  is.cal <- is(origpos, "positionsCalendar")
  ## construct an object for return
  if(is.cal) {
    newobj <- timeSeries()
    newobj@fiscal.year.start <- x@fiscal.year.start
    newobj@type <- x@type
  }
  else {
    newobj <- signalSeries()
    newobj@units.position <- x@units.position
  }
  newobj@units <- x@units
  newobj@title <- x@title
  newobj@documentation <- x@documentation
  newobj@attributes <- x@attributes
  if(length(origpos) < 1)
    return(newobj)
  rng <- range(origpos)
  if(moving) {
    if(!missing(pos))
      warning("pos argument ignored for moving aggregation")
    len <- length(origpos)
    if(len < moving)
      return(newobj)
    ## set up positions for returns and offsets
    poslen <- trunc((len - moving)/incr) + 1
    startindx <- (0:(poslen - 1)) * incr + 1
    pos <- origpos[startindx]
    endindx <- startindx + moving - 1
    diffpos <- as(origpos, "numeric")
    diffpos <- diffpos[endindx] - diffpos[startindx]
    whichpos <- TRUE
    ## call aggregation function
    newdat <- lapply(startindx,
                     function(rowstart, numrows, x, nc, func, together, ...)
                     {
                       i = rowstart:(rowstart +  numrows - 1)
                       cols <- if(together){
                         func(subscript2d(x,i,))
                       } else if(nc > 0)
                         sapply( 1:nc, function(col, x, fun, ...){
                           fun(subscript2d(x,,col)[,1] )
                         },
                                subscript2d(x,i,),
                                func)
                       else numeric(0)
                       cols <- as.rectangular(cols)
                       if(is(cols, "matrix") && numRows(cols) != 1)
                         cols <- matrix(cols, nrow = 1)
                       cols
                     }
                     , moving, x@data, ncol(x), FUN, together, ...)
    newdat <- do.call("rbind", newdat)
  } else {
    ## construct positions if missing
    if(missing(pos) && is.cal) {
      ## make sure to go just beyond top end of range
      if(is(by, "character")) {
        endpt <- rng[2]
        endpt@columns[[2]] <- endpt@columns[[2]] + 1L
        pos <- timeSeq(from = rng[1], to = endpt, by = 
                       by, k.by = k.by, align.by = align.by,
                       extend = TRUE, week.align = week.align,
                       holidays = holidays)
      }
      else pos <- timeSeq(from = rng[1], to = rng[2] + by,
                          by = by)
      diffpos <- diff(as(pos, "numeric"))
    }
    else if(missing(pos)) {
      pos <- seq(from = rng[1], to = rng[2] + by, by = by)
      diffpos <- diff(as(pos, "numeric"))
    }
    else {
      ## if they passed in positions, add another at the end,
      ## to make into breakpoints
      poslen <- length(pos)
      diffpos <- diff(as(pos, "numeric"))
      if(any(is.na(diffpos)) || min(diffpos) < 0)
        stop("Positions for aggregation must be monotonically increasing without NA values"
             )
      if(poslen > 1) {
        dm <- max(diffpos)
        pos[poslen + 1] <- pos[poslen] + dm
        diffpos <- c(diffpos, dm)
      }
      else {
        pos[poslen + 1] <- pos[poslen] + 1
        diffpos <- c(diffpos, 1)
      }
    }
    poslen <- length(pos)
    if(length(pos) < 2)
      return(newobj)
    ## make bins, and take endpoint off positions we are keeping
    posbins <- pos
    pos <- pos[ - length(pos)]
    ## call cut to figure out which positions go in which bins
    ## figure out the next diff in sequence, if possible
    ## get right bin endpoints
    if(include.ends) {
      posbins[1] <- min(posbins[1], rng[1])
      posbins[poslen] <- max(posbins[poslen], rng[2])
    }
    poscut <- cut(origpos, posbins, include.lowest = FALSE, right=FALSE)
    ## poscut[is.na(poscut)] <- 0
    if(!drop.empty) {
      poscats <- seq(along = pos)
      whichpos <- TRUE
    }
    else {
      poscats <- sort(unique(poscut))
      ## poscats <- poscats[poscats > 0]
      whichpos <- poscats
    }
    if(is.character(FUN) && substring(FUN, 1, 4) == "fast") {
      ## fast aggregation
      FUN <- get(paste("igroup", substring(FUN, 5), sep = "")
                 )
      nc <- numCols(x)
      if(nc == 1)
        newdat <- FUN(x@data, poscut, ...)
      else {
        newdat <- lapply(1:nc, function(col, dat, grps, fun, ...)
                         {
                           fun(dat[,col], grps, ...)
                         }
                         , x@data, poscut, FUN, ...)
        newdat <- do.call("data.frame", newdat)
      }
    }
    else {
      ## general aggregation
      ## apply the function to get new series data
      newdat <- lapply(poscats, function(thecat, x, nc, 
                                         poscut, func, together, ...)
                       {
                         cols <- if(together)
                           func(subscript2d(x, poscut == thecat,))
                         else if(nc > 0)
                           sapply(1:nc, function(col, x, fun, ...)
                                  {
                                    if( !is.null(dim(x)) )
                                      x <- subscript2d(x,, col)[,1]
                                    fun( x[!is.na(x)] )
                                  }
                                  , subscript2d(x, poscut == thecat,),
                                  func, ...)
                         else numeric(0)
                         cols <- as.rectangular(cols)
                         if(is(cols, "matrix") && numRows(cols) != 1)
                           cols <- matrix(cols, nrow = 1)
                         cols
                       }
                       , x@data, ncol(x), poscut, FUN, together, ...)
      newdat <- do.call("rbind", newdat)
    }
  }
  ## add offset or adj to positions
  if(!missing(offset)) pos <- pos + offset
  else if(!missing(adj) && (poslen > 1)) {
		pos <- pos + adj * diffpos
              }
  if(missing(colnames))
    colnames <- names(x@data)
  if(!is.null(colnames) && (numCols(newdat) == length(colnames)))
    colnames(newdat) <- colnames
  newobj@positions <- pos[whichpos]
  newobj@data <- newdat
  newobj
}

"hloc" <- 
function(x)
{
  ## return the high, low, open, and closing values of a vector in that order
  if(!is.atomic(x)) {
      stop("x must be atomic")
  }
  if(length(dim(x))) x <- as.vector(as.matrix(x))
  len <- length(x)
  if(!len)
    as(rep(NA, 4), class(x))
  else c(range(x)[c(2, 1)], x[c(1, len)])
}


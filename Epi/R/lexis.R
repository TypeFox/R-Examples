Lexis <-
function(entry, exit, duration, entry.status=0, exit.status=0, id, data,
         merge=TRUE, states, tol=.Machine$double.eps^0.5,
         keep.dropped=FALSE )
{
  nmissing <- missing(entry) + missing(exit) + missing(duration)
  if (nmissing > 2)
        stop("At least one of the arguments exit and duration must be supplied")
  only.exit <- missing( entry.status ) && !missing( exit.status )

  ## If data argument is supplied, use it to evaluate arguments
  if (!missing(data)) {
    if (!missing(entry)) {
      entry <- eval(substitute(entry), data, parent.frame())
    }
    if (!missing(exit)) {
      exit <- eval(substitute(exit), data, parent.frame())
    }
    if (!missing(duration)) {
      duration <- eval(substitute(duration), data, parent.frame())
    }
    entry.status <- eval(substitute(entry.status), data, parent.frame())
    exit.status <- eval(substitute(exit.status), data, parent.frame())
    if (!missing(id)) {
      id <- eval(substitute(id), data, parent.frame())
    }
    if (merge) {
      data <- as.data.frame(data)
    }
  }

  ## Check for missing values in status variables
  wh.miss <- any(is.na(entry.status)) + 2*any(is.na(exit.status))
  if ( wh.miss > 0 ) stop("Missing values in ",
                          switch( wh.miss, "entry status",
                                           "exit status",
                                           "entry AND exit status" ) )

  ## Adjust entry status mode according to exit status
  if( only.exit )
    {
    if( is.logical( exit.status ) )
        entry.status <- FALSE
    if( is.character( exit.status ) )
        exit.status <- factor( exit.status )
    if( is.factor( exit.status ) )
        {
        entry.status <- factor( rep( levels(exit.status)[1],
                                     length(exit.status)),
                                levels=levels(exit.status),
                                labels=levels(exit.status) )
        cat("NOTE: entry.status has been set to",
            paste( '"', levels(exit.status)[1], '"', sep='' ),
            "for all.\n" )
        }
    if( is.numeric( exit.status ) )
        entry.status <- rep( 0, length( exit.status ) )
    }

  ## Convert character states to factors
  if( is.character(entry.status) ) entry.status <- factor(entry.status)
  if( is.character( exit.status) )  exit.status <- factor( exit.status)

  ## Check compatibility of entry and exit status
  if (is.factor(entry.status) || is.factor(exit.status)) {
      if (is.factor(entry.status) && is.factor(exit.status)) {
          if (!identical(levels(entry.status),levels(exit.status))) {
              all.levels = union(levels(entry.status),levels(exit.status))
              entry.status <- factor( entry.status, levels=all.levels )
               exit.status <- factor(  exit.status, levels=all.levels )
            cat("Incompatible factor levels in entry.status and exit.status:\n",
                "both lex.Cst and lex.Xst now have levels:\n", all.levels, "\n")
          }
      }
      else {
          stop("Incompatible classes for entry and exit status")
      }
  }
  else {
      if (mode(entry.status) != mode(exit.status)) {
          stop("Incompatible mode for entry and exit status")
      }
  }

  ## If entry is missing and one of the others is given as a list of length
  ## one, entry is assumed to be 0 on this only timescale.
  if( nmissing==2 )
    {
    if( !missing(exit) )
      { if( length(exit)>1 )
          stop("If 'entry' is omitted, only one timescale can be specified.")
        else
        {
        entry <- exit
        entry[[1]] <- 0*entry[[1]]
        cat( "NOTE: entry is assumed to be 0 on the",names(exit),"timescale.\n")
        }
      }
    else
    if( !missing(duration) )
      { if( length(duration)>1 )
          stop("If 'entry' is omitted, only one timescale can be specified")
        else
        {
        entry <- duration
        entry[[1]] <- 0*entry[[1]]
        cat( "NOTE: entry is assumed to be 0 on the",names(duration),"timescale.\n")
        }
      }
    else
    stop("Either exit or duration must be supplied.")
    }

  ## Coerce entry and exit lists to data frames

  if(!missing(entry)) {
    entry <- as.data.frame(entry)
    if (is.null(names(entry)))
      stop("entry times have no names")
    if (any(substr(names(entry),1,4) == "lex."))
      stop("names starting with \"lex.\" cannot be used for time scales")
  }

  if(!missing(exit)) {
    exit <- as.data.frame(exit)
    if (is.null(names(exit)))
      stop("exit times have no names")
    if (any(substr(names(exit),1,4) == "lex."))
      stop("names starting with \"lex.\" cannot be used for time scales")
  }

  if(!missing(duration)) {
    duration <- as.data.frame(duration)
    if (is.null(names(duration)))
      stop("duration have no names")
    if (any(substr(names(duration),1,4) == "lex."))
      stop("names starting with \"lex.\" cannot be used for time scales")
  }

  if (missing(entry)) {
    ## Impute entry
    entry <- exit - duration
  }

  if (missing(duration)) {
    ## Impute duration
    full.time.scales <- intersect(names(entry), names(exit))
    if (length(full.time.scales) == 0) {
      stop("Cannot calculate duration from entry and exit times")
    }
    duration <- exit[,full.time.scales[1]] - entry[,full.time.scales[1]]
  }

  if (missing(exit)) {
    all.time.scales <- names(entry)
  }
  else {
    ## We dont need the exit times but, if they are supplied, we must
    ## make sure they are consistent with the entry and duration.

    all.time.scales <- unique(c(names(entry), names(exit)))
    ## Fill in any missing entry times
    entry.missing <- setdiff(all.time.scales, names(entry))
    if (length(entry.missing) > 0) {
      entry <- cbind(entry, exit[,entry.missing, drop=FALSE] - duration)
      }
    ## Check that duration is the same on all time scales
    dura <- exit - entry[,names(exit),drop=FALSE]
    if (missing(duration)) {
      duration <- dura[,1] #BxC# apply( dura, 1, mean, na.rm=TRUE )
                           # Allows for timescales with missing values
      }
    ok <- sapply(lapply(dura, all.equal, duration), isTRUE)
#   ok <- sapply(lapply(dura, all.equal, duration),
#                function(x) identical(FALSE,x) )
    if (!all(ok)) {
      stop("Duration is not the same on all time scales")
    }
  }
#  Taken care of by the code that detects whether lex.du <= tol
#  ## Check that duration is positive
#  if (any(duration<0)) {
#    stop("Duration must be non-negative")
#  }

  ## Make sure id value - if supplied - is valid. Otherwise supply default id

  if (missing(id)) {
    id <- 1:nrow(entry)
  }
  else if (any(duplicated(id))) {
    ##Fixme: check for overlapping intervals
    ##stop("Duplicate values in id")
  }

  ## Return a data frame with the entry times, duration, and status
  ## variables Use the prefix "lex." for the names of reserved
  ## variables.
  if( is.data.frame( duration ) ) duration <- duration[,1]
  lex <- data.frame(entry, "lex.dur" = duration,
                           "lex.Cst" = entry.status,
                           "lex.Xst" = exit.status,
                           "lex.id"  = id )

  #### Addition by BxC --- support for states as factors
  # Convert states to factors if states are given
  if( !missing( states ) ) #is.character( states ) )
    {
    # This as.character-business is necessary because we cannot assume
    # that the values of states are 1,2, etc.
    st.lev <- sort( unique( as.character( c(lex$lex.Cst,lex$lex.Xst) ) ) )
    lex$lex.Cst <- factor( as.character(lex$lex.Cst), levels=st.lev, labels=states )
    lex$lex.Xst <- factor( as.character(lex$lex.Xst), levels=st.lev, labels=states )
    }

  if (!missing(data) && merge) {
    duplicate.names <- intersect(names(lex), names(data))
    if (length(duplicate.names) > 0) {
      stop("Cannot merge data with duplicate names")
    }
    lex <- cbind(lex, data)
  }

  ## Drop rows with short or negantive duration for consistency with splitLexis
  short.dur <- lex$lex.dur <= tol
  if ( any(short.dur) ) {
      warning("Dropping ", sum(short.dur),
              " rows with duration of follow up < tol\n",
      if( keep.dropped ) "  The dropped rows are in the attribute 'dropped'\n",
      if( keep.dropped ) "  To see them type attr(Obj,'dropped'),\n",
      if( keep.dropped ) "  to get rid of them type attr(Obj,'dropped') <- NULL\n",
      if( keep.dropped ) "  - where 'Obj' is the name of your Lexis object" )
      lex <- subset(lex, !short.dur)
      if( keep.dropped ) attr(lex,"dropped") <- subset(data, short.dur)
      }

  ## Return Lexis object
  attr(lex,"time.scales") <- all.time.scales
  attr(lex,"time.since") <- rep( "", length(all.time.scales) )
  breaks <- vector("list", length(all.time.scales))
  names(breaks) <- all.time.scales
  attr(lex,"breaks") <- breaks
  class(lex) <- c("Lexis", class(lex))
  return(lex)
}

is.Lexis <- function(x)
{
  inherits(x, "Lexis")
}

check.time.scale <- function(lex, time.scale=NULL)
{

  ## Utility function, returns the names of the time scales in a Lexis object
  ## lex - a Lexis object
  ## time.scale - a numeric or character vector. The function checks that
  ##              these are valid time scales for the Lexis object.
  ## Return value is a character vector containing the  names of the requested
  ## time scales

  all.names <- timeScales(lex)
  if (is.null(time.scale))
    return(all.names)

  nscale <- length(time.scale)
  scale.names <- character(nscale)
  if (is.character(time.scale)) {
    for (i in 1:nscale) {
      if (is.null(lex[[time.scale[i]]]))
        stop(time.scale[i], " is not a valid time scale name")
    }
  }
  else if (is.numeric(time.scale)) {
    if (any(time.scale > length(all.names))  || any(time.scale < 1))
      stop(time.scale, " not valid time scale column number(s)")
    time.scale <- all.names[time.scale]
  }
  else {
    stop("invalid type for time scale")
  }
  return(time.scale)
}

valid.times <-
function(x, time.scale=1)
  {
  # A utility function that returns a data.frame / Lexis object with
  # rows with missing timescales removed
  x[complete.cases(x[,check.time.scale(x,time.scale)]),]
  }

plot.Lexis.1D <- function(x, time.scale=1, breaks="lightgray",
                          type="l", col="darkgray", xlim, ylim, xlab, ylab,
                          ...)
{
  ## x Lexis object
  ## time.scale  name of time scale to plot

  if (length(time.scale) != 1)
    stop("Only one time scale allowed")

  x <- valid.times(x,time.scale)
  time.entry <- x[,time.scale]
  time.exit <- x[,time.scale] + x$lex.dur
  id <- x$lex.id

  if (missing(xlim))
    xlim <- c(min(time.entry), max(time.exit))
  if (missing(ylim))
    ylim <- range(id)
  if (missing(xlab))
    xlab <- time.scale
  if (missing(ylab))
    ylab <- "id number"

  plot(time.entry, id, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
       ...)
  if (type=="b" || type=="l") {
    segments(time.entry, id, time.exit, id, col=col, ...)
  }
  if (type=="b" || type=="p") {
    points(time.exit, id, col=col, ...)
  }
  ## Plot break points
  brk <- attr(x,"breaks")[[time.scale]]
  abline(v=brk, col=breaks, ...)
}

points.Lexis.1D <- function(x, time.scale, ...)
{
  x <- valid.times(x,time.scale)
  time.exit <- x[,time.scale] + x$lex.dur
  points(time.exit, x$lex.id, ...)
}

lines.Lexis.1D <- function(x, time.scale, type="l", col="darkgray", breaks="lightgray",
                           ...)
{
  x <- valid.times(x,time.scale)
  time.entry <- x[,time.scale]
  time.exit <- x[,time.scale] + x$lex.dur
  id <- x$lex.id
  segments(time.entry, id, time.exit, id, col=col, ...)
  ## Plot break points
  brk <- attr(x,"breaks")[[time.scale]]
  abline(v=brk, col=breaks, ...)

}

plot.Lexis.2D <- function(x, time.scale, breaks="lightgray",
                          type="l", col="darkgray",
                          xlim, ylim, xlab, ylab,
                          grid=FALSE,
                      col.grid="lightgray",
                      lty.grid=2,
                      coh.grid=FALSE,
                          ...)
{
  if (length(time.scale) != 2)
    stop("Two time scales are required")

  x <- valid.times(x,time.scale)

  time.entry <- time.exit <- vector("list",2)
  for (i in 1:2) {
    time.entry[[i]] <- x[,time.scale[i]]
    time.exit[[i]] <- x[,time.scale[i]] + x$lex.dur
  }

  if (missing(xlim) && missing(ylim)) {
    ## If no axis limits are given, set the plotting region to be
    ## square, and adjust the axis limits to cover the same time interval.
    ## All life lines will then be at 45 degrees.
    opar <- par(pty="s")
    on.exit(par(opar))
    min.times <- sapply(time.entry, min)
    max.times <- sapply(time.exit, max)
    xywidth <- max(max.times - min.times)
    xlim <- min.times[1] + c(0, xywidth)
    ylim <- min.times[2] + c(0, xywidth)
  }
  else if (missing(xlim)) {
    xlim <- c(min(time.entry[[1]]), max(time.exit[[1]]))
  }
  else if (missing(ylim)) {
    ylim <- c(min(time.entry[[2]]), max(time.exit[[2]]))
  }

  if (missing(xlab))
    xlab <- time.scale[1]
  if (missing(ylab))
    ylab <- time.scale[2]

  plot(time.entry[[1]], time.entry[[2]], type="n",
       xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)

# Set up the background grid(s):
  if (!missing(grid)) {
     if (is.logical(grid)) {
        if (grid) {
        vgrid <- pretty(xlim)
        hgrid <- pretty(ylim)
        } }
     else if (is.list(grid)) {
        vgrid <- grid[[1]]
        hgrid <- grid[[length(grid)]]
        }
     else if (is.numeric(grid)) {
       vgrid <- grid - min( grid ) + min( pretty( xlim )[pretty(xlim)>=par("usr")[1]] )
       hgrid <- grid - min( grid ) + min( pretty( ylim )[pretty(ylim)>=par("usr")[3]] )
        }
     else stop( "'grid' must be either logical, list or a numeric vector" )
   # and plot the grid:
     abline( v=vgrid, h=hgrid, col=col.grid, lty=lty.grid )
     box()
     }
  if (!missing(grid) & coh.grid) {
     # Make the 45-degree grids as fine as the finest grid on the axes
     for (yy in c(hgrid-diff(range(hgrid)),hgrid))
         abline( yy-min(vgrid), 1, col=col.grid, lty=lty.grid )
     for (yy in c(vgrid-diff(range(vgrid)),vgrid))
         abline( min(hgrid)-yy, 1, col=col.grid, lty=lty.grid )
     }
# End of explicitly requested background grid(s) (PHEW!)

  if (type=="b" || type=="l") {
    segments(time.entry[[1]], time.entry[[2]], time.exit[[1]], time.exit[[2]],
             col=col, ...)
    }
  if (type=="b" || type=="p") {
    points(time.exit[[1]], time.exit[[2]], col = col, ...)
  }
  if (type != "n") {
    ## Plot break points
    brk <- attr(x,"breaks")[time.scale]
    abline(v=brk[[1]], h=brk[[2]], col=breaks, ...)
  }
}

points.Lexis.2D <- function(x, time.scale, ...)
{
  x <- valid.times(x,time.scale)
  time.exit <- vector("list",2)
  for (i in 1:2) {
    time.exit[[i]] <- x[,time.scale[i]] + x$lex.dur
  }
  points( time.exit[[1]], time.exit[[2]], ...)
}

lines.Lexis.2D <- function(x, time.scale, col="darkgray", ...)
{
  x <- valid.times(x,time.scale)
  time.entry <- time.exit <- vector("list",2)
  for (i in 1:2) {
    time.entry[[i]] <- x[,time.scale[i]]
    time.exit[[i]] <- x[,time.scale[i]] + x$lex.dur
  }
  segments(time.entry[[1]], time.entry[[2]], time.exit[[1]], time.exit[[2]],
           col=col, ...)
}

### Plotting generic functions
plot.Lexis <-
function( x = Lexis( entry=list(Date=1900,Age=0), exit=list(Age=0) ),
          time.scale=NULL, type="l", breaks="lightgray", ...)
{
  time.scale <- check.time.scale(x, time.scale)
  if (length(time.scale) > 2)
    time.scale <- time.scale[1:2]
  # Save the timescale(s) for use in subsequent calls
  options( Lexis.time.scale = time.scale )

  if (length(time.scale) == 1)
    plot.Lexis.1D(x, time.scale=time.scale, type=type, breaks=breaks, ...)
  else if (length(time.scale) == 2)
    plot.Lexis.2D(x, time.scale=time.scale, type=type, breaks=breaks, ...)
}


lines.Lexis <- function(x, time.scale=options()[["Lexis.time.scale"]], ...)
{
  time.scale <- check.time.scale(x, time.scale)
  if (length(time.scale) > 2)
    time.scale <- time.scale[1:2]

  if (length(time.scale) == 1)
    lines.Lexis.1D(x, time.scale=time.scale, ...)
  else if (length(time.scale) == 2)
    lines.Lexis.2D(x, time.scale=time.scale, ...)
}

points.Lexis <- function(x, time.scale=options()[["Lexis.time.scale"]], ...)
{
  time.scale <- check.time.scale(x, time.scale)
  if (length(time.scale) > 2)
    time.scale <- time.scale[1:2]
  if (length(time.scale) == 1)
    points.Lexis.1D(x, time.scale=time.scale, ...)
  else if (length(time.scale) == 2)
    points.Lexis.2D(x, time.scale=time.scale, ...)
}

PY.ann <- function (x, ...) UseMethod("PY.ann")
PY.ann.Lexis <-
function( x, time.scale=options()[["Lexis.time.scale"]], digits=1, ... )
{
  if( !inherits(x,"Lexis") )
    stop( "Only meaningful for Lexis objects not for objects of class ", class(x) )
  wh.x <- x[,time.scale[1]] + x[,"lex.dur"]/2
  if( two.scales <- length(time.scale)==2 )
     wh.y <- x[,time.scale[2]] + x[,"lex.dur"]/2
  else
     wh.y <- x[,"lex.id"]
  text( wh.x, wh.y,
        formatC(x$lex.dur,format="f",digits=digits),
        adj=c(0.5,-0.5), srt=if(two.scales) 45 else 0, ... )
}

### Generic functions

### Methods for data.frame drop Lexis attributes, so we need a Lexis
### method that adds them again

subset.Lexis <- function(x, ...)
{
  y <-  subset.data.frame(x, ...)
  attr(y,"breaks") <- attr(x, "breaks")
  attr(y,"time.scales") <- attr(x, "time.scales")
  attr(y,"time.since") <- attr(x, "time.since")
  return(y)
}

`[.Lexis` <-
function( x, ... )
{
    structure( NextMethod(),
              breaks = attr(x, "breaks"),
              time.scales = attr(x, "time.scales"),
              time.since = attr(x, "time.since") )
}

merge.data.frame <- function(x, y, ...)
{
  if (is.Lexis(x))
    merge.Lexis(x, y, ...)
  else if (is.Lexis(y))
    merge.Lexis(y, x, ...)
  else
    base::merge.data.frame(x, y, ...)
}

merge.Lexis <- function(x, y, id, by, ...)
{
  if (!missing(id)) {
    if (!is.character(id) || length(id) != 1 || !(id %in% names(y))) {
      stop("id must be the name of a single variable in y")
    }
    if (any(duplicated(y[[id]]))) {
      stop("values of the id variable must be unique in y")
    }
    y$lex.id <- y[[id]]
  }
  else if (missing(by)) {
    by <- intersect(names(x), names(y))
    if (length(by)==0) {
      stop("x and y have no variable names in common")
    }
  }

  z <-  base::merge.data.frame(x, y, ...)
  attr(z,"breaks") <- attr(x, "breaks")
  attr(z,"time.scales") <- attr(x, "time.scales")
  attr(z,"time.since") <- attr(x, "time.since")
  class(z) <- c("Lexis", "data.frame")
  return(z)
}

cbind.Lexis <-
function( ... )
{
allargs <- list( ... )
# Check that at least one argument is Lexis
is.lex <- sapply( allargs, inherits, "Lexis" )
if( all(!is.lex) ) stop( "At least one argument nust be a Lexis object\n",
                         "and none of the given are.\n")
if( sum(is.lex)>1 ) stop( "It is meaningless to 'cbind' several Lexis objects:",
                          " arguments ", paste( which(is.lex), collapse=","),
                          " are Lexis objects.\n" )
is.lex <- which(is.lex) 
res <- do.call( base::cbind.data.frame, allargs )
attr( res, "class"       ) <- attr( allargs[[is.lex]], "class"       )
attr( res, "breaks"      ) <- attr( allargs[[is.lex]], "breaks"      )
attr( res, "time.scales" ) <- attr( allargs[[is.lex]], "time.scales" )
attr( res, "time.since"  ) <- attr( allargs[[is.lex]], "time.since"  )
res
}

rbind.Lexis <-
function( ... )
{
# A list of all Lexis objects    
allargs <- list( ... )
# Check if they are all Lexis
# (or possibly NULL - often rbind-ing with NULL is very useful)
is.lex <- sapply( allargs, inherits, "Lexis" )
is.nul <- sapply( allargs, is.null )
if( !all(is.lex[!is.nul]) )
    stop( "All arguments must be Lexis objects,\n",
          "arguments number ", which(!is.lex & !is.nul), " are not." )
# Put them all together
allargs <- allargs[!is.nul]
res <- plyr::rbind.fill( allargs )
# Get the union of time.scale names and the corresponding time.since
tscl <- do.call( c, lapply( allargs, function(x) attr(x,"time.scales") ) )
tsin <- do.call( c, lapply( allargs, function(x) attr(x,"time.since" ) ) )
# but only one copy of each
scls <- match( unique(tscl), tscl )
tscl <- tscl[scls]
tsin <- tsin[scls]    
# Fish out the breaks on timescale in turn from all input objects and
# - if all the non-NULL are identical use this
# - if not, set the corresponding break to NULL
newbrks <- list()
# all the breaks attributes in a list    
brks <- lapply( allargs, function(x) attr(x,"breaks") )
# run through the timescales found 
for( scl in tscl )
   {
   # breaks for this timescale in any of the objects   
   brk <- lapply( brks, function(x) x[[scl]] )
   # but only the non-null ones
   brk <- brk[!sapply( brk, is.null )]
   # if more than one occurrence, all non-NULL breaks should be identical
   if( ( length(brk)>1 & 
         all( sapply( brk[-1], function(x) identical(brk[[1]],x) ) ) )
       | length(brk) == 1 ) newbrks[scl] <- brk[1]
   else newbrks[scl] <- list(NULL)    
   }
# define attributes of the reulting object:
attr( res, "class"       ) <- c( "Lexis", "data.frame" )
attr( res, "breaks"      ) <- newbrks
attr( res, "time.scales" ) <- tscl
attr( res, "time.since"  ) <- tsin
res  
}

## Extractor functions

entry <- function(x, time.scale = NULL, by.id = FALSE )
{
    time.scale <- check.time.scale(x, time.scale)
    wh <- x[,time.scale[1]] ==
     ave( x[,time.scale[1]], x$lex.id, FUN=if( by.id ) min else I )
    if (length(time.scale) > 1) {
        res <- as.matrix(x[wh, time.scale])
        if( by.id ) rownames( res ) <- x$lex.id[wh]
        return( res )
    }
    else {
        res <- x[wh, time.scale]
        if( by.id ) names( res ) <- x$lex.id[wh]
        return( res )
    }
}

exit <- function(x, time.scale = NULL, by.id = FALSE )
{
    time.scale <- check.time.scale(x, time.scale)
    wh <- x[,time.scale[1]] ==
     ave( x[,time.scale[1]], x$lex.id, FUN=if( by.id ) max else I )
    if (length(time.scale) > 1) {
        res <- as.matrix(x[wh, time.scale]) + x$lex.dur[wh]
        if( by.id ) rownames( res ) <- x$lex.id[wh]
        return( res )
    }
    else {
        res <- x[wh, time.scale] + x$lex.dur[wh]
        if( by.id ) names( res ) <- x$lex.id[wh]
        return( res )
    }
}

dur <- function(x, by.id=FALSE)
{
  if( by.id ) return( tapply(x$lex.dur,x$lex.id,sum) )
  else        return(        x$lex.dur               )
}

status <- function(x, at="exit", by.id = FALSE)
{
  at <- match.arg(at, c("entry","exit"))
  wh <- x[,timeScales(x)[1]] ==
   ave( x[,timeScales(x)[1]], x$lex.id, FUN=if(by.id)
                                             switch(at,
                                                    "entry"=min,
                                                     "exit"=max)
                                             else I )
  res <- switch(at, "entry"=x$lex.Cst, "exit"=x$lex.Xst)[wh]
  if( by.id ) names( res ) <- x$lex.id[wh]
  res
}

timeScales <- function(x)
{
  return (attr(x,"time.scales"))
}

timeBand <- function(lex, time.scale, type="integer")
{
  time.scale <- check.time.scale(lex, time.scale)[1]
  breaks <- attr(lex, "breaks")[[time.scale]]

  time1 <- lex[[time.scale]]
  band <- findInterval(time1, breaks)

  ##Check that right hand side of interval falls in the same band
  abrk <- c(breaks, Inf)
  tol <- sqrt(.Machine$double.eps)
  if (any(time1 + lex$lex.dur > abrk[band+1] + tol)) {
    stop("Intervals spanning multiple time bands in Lexis object")
  }

  type <- match.arg(type, choices = c("integer","factor","left","middle",
                    "right"))
  if (type=="integer") {
     return(band)
  }

  I1 <- c(-Inf, breaks)
  I2 <- c(breaks, Inf)
  labels <- switch(type,
                   "factor" = paste("(", I1, ",", I2, "]", sep=""),
                   "left" = I1,
                   "right" = I2,
                   "middle" = (I1 + I2)/2)

  if(type=="factor") {
    return(factor(band, levels=0:length(breaks), labels=labels))
  }
  else {
    return(labels[band+1])
  }
}

breaks <- function(lex, time.scale)
{
  time.scale <- check.time.scale(lex, time.scale)[1]
  return(attr(lex, "breaks")[[time.scale]])
}

transform.Lexis <- function(`_data`, ... )
{
    save.at <- attributes(`_data`)
    ## We can't use NextMethod here because of the special scoping rules
    ## used by transform.data.frame
    y <- base::transform.data.frame(`_data`, ...)
    save.at[["names"]] <- attr(y, "names")
    attributes(y) <- save.at
    y
}

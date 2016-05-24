is.tframed.its <- function(x) {TRUE}

tframe.its <- function (x) {
  tf <- its::dates(x)
  class(tf) <- c( "itstframe", class(tf), "tframe")
  tf
  }

tfUnSet.its <- function(x)      {x@.Data}
tfSet.itstframe <- function(value, x) {
   class(value) <- class(value)[class(value) != "itstframe"]
   its::its(x, value)
   }

"seriesNames<-.its" <- function (x, value) 
  {if (is.matrix(x)) dimnames(x) <- list(NULL, value)
   else attr(x, "seriesNames") <- value
   x
  }

tfstart.its <- function(x) its::dates(x)[1] # start(x) returns character rather than a date
tfend.its   <- function(x) its::dates(x)[Tobs(x)]
time.its    <- function(x, ...) its::dates(x, ...)

Tobs.its <- function(x)  NROW(x)

tfstart.itstframe <- function(x) x[1]
tfend.itstframe   <- function(x) x[length(x)]
Tobs.itstframe     <- function(x) length(x)

tfwindow.its <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
  {# With the default warn=T warnings will be issued if no truncation takes
   #  place because start or end is outside the range of data.
   if (!warn) 
     {opts <- options(warn = -1)
      on.exit(options(opts))
     }
   d <- its::dates(x)
   i <- rep(TRUE, Tobs(x))
   if(!is.null(start)) i <- i &  
             (d >= if(is.character(start)) as.POSIXct(start) else start)
   if(!is.null(end))   i <- i &  
             (d <= if(is.character(end)) as.POSIXct(end) else end)
   y <- its::its(x[i,, drop=FALSE], d[i])
   seriesNames(y) <- seriesNames(x)
   y
  }

tbind.its <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
 {nm <- seriesNames(x)
  for (z in list(...)) {
    if (!is.null(z)) {
      nm <- c(nm, seriesNames(z))
      x <- its::union(x, z)
      }
    }
  if (!pad.start | !pad.end)
     x <- trimNA(x, startNAs= !pad.start, endNAs= !pad.end)
  seriesNames(x) <- nm
  x
 }  


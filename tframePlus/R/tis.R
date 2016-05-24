is.tframed.tis <- function(x) {TRUE}

tframe.tis <- function (x) {
  tf <- time(x)
  attr(tf, "start") <- attr(x, "start")
  class(tf) <- c( "tistframe", class(tf), "tframe")
  tf
  }

tfSet.tistframe <- function(value, x) {
   r <- tis::tis(x, start=attr(value, "start") )
   if (inherits(r, "try-error")) {r <- x ; attr(r, "tframe") <- value}
   r
   }

tfstart.tis <- function(x) start(x)  #tis::start(x)
tfend.tis   <- function(x) end(x)    #tis::end(x)

Tobs.tis <- function(x)  NROW(x)

tfstart.tistframe <- function(x) x[1]
tfend.tistframe   <- function(x) x[length(x)]
Tobs.tistframe     <- function(x) length(x)

tfwindow.tis <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
  {# With the default warn=T warnings will be issued if no truncation takes
   #  place because start or end is outside the range of data.
   y <- window(x, start=start, end=end, noWarn=!warn)
   seriesNames(y) <- seriesNames(x)
   y
  }

tbind.tis <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
 {nm <- seriesNames(x)
  for (z in list(...)) {
    if (!is.null(z)) {
      nm <- c(nm, seriesNames(z))
      x <- cbind(x, z)
      }
    }
  if (!pad.start | !pad.end)
     x <- trimNA(x, startNAs= !pad.start, endNAs= !pad.end)
  seriesNames(x) <- nm
  x
 }  


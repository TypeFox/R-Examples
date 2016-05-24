is.tframed.xts <- function(x) {TRUE}

tframe.xts <- function (x) {
  tf <- zoo::index(x)
  class(tf) <- c( "xtstframe", class(tf), "tframe")
  tf
  }

tfUnSet.xts <- function(x)  zoo::coredata(x)

tfSet.xtstframe <- function(value, x) { 
  class(value) <- class(value)[class(value) != "xtstframe"]
  r <- xts::xts(x, value)  
  seriesNames(r) <- seriesNames(x)
  r
 }

"seriesNames<-.xts" <- function (x, value) 
  {if (is.matrix(x)) dimnames(x) <- list(NULL, value)
   else attr(x, "seriesNames") <- value
   x
  }

Tobs.xts <- function(x)  NROW(x)

tfstart.xtstframe <- function(x) x[1]
tfend.xtstframe   <- function(x) x[length(x)]
Tobs.xtstframe     <- function(x) length(x)

tfL.xts <- function (x, p = 1) lag(x, k = -p)

tfwindow.xts <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
  {# With the default warn=T warnings will be issued if no truncation takes
   #  place because start or end is outside the range of data.
   if (!warn) 
     {opts <- options(warn = -1)
      on.exit(options(opts))
     }
   y <- window(x, start=start, end=end)
   seriesNames(y) <- seriesNames(x)
   attr(y, "TSrefperiod") <- attr(x, "TSrefperiod")
   y
  }

tfExpand.xts <- function(x, add.start = 0, add.end = 0){
   idx <- time(x)
   r <- as.matrix(zoo::coredata(x))
   if (add.start > 0 ) {
     idx <- c(start(x) - seq(add.start), idx)
     r <- rbind(matrix(NA, add.start, ncol(r)), r)
     }
   if (add.end > 0 ) {
     idx <- c(idx, end(x) + seq(add.end))
     r <- rbind(r, matrix(NA,add.end, ncol(r)))
     }
   xts::xts(r, order.by = idx) 
   }

tbind.xts <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
 {nm <- seriesNames(x)
  ref <- attr(x, "TSrefperiod")
  for (z in list(...)) {
    if (!is.null(z)) {
      nm  <- c(nm,  seriesNames(z))
      ref <- c(ref, attr(z, "TSrefperiod"))
      x <- cbind(x, z)
      }
    }
  if (!pad.start | !pad.end)
     x <- trimNA(x, startNAs= !pad.start, endNAs= !pad.end)
  seriesNames(x) <- nm
  attr(x, "TSrefperiod") <- ref
  x
 }  


is.tframed.zoo <- function(x) {TRUE}

tframe.zoo <- function (x) {
  tf <- zoo::index(x)
  class(tf) <- c( "zootframe", class(tf), "tframe")
  tf
  }

tfUnSet.zoo <- function(x)      {zoo::coredata(x)}
tfSet.zootframe <- function(value, x){ 
  if(Tobs(value) != Tobs(x)) stop("number of Tobs of observations must correspond to number of Tobs indicated by tframe.")
  class(value) <- class(value)[class(value) != "zootframe"]
  r <- zoo::zoo(x, order.by = value) 
  seriesNames(r) <- seriesNames(x)
  r
  }

tfSet.Date    <- function(value, x){ 
  if(Tobs(value) != Tobs(x)) stop("number of Tobs of observations must correspond to number of Tobs indicated by tframe.")
  r <- zoo::zoo(x, order.by = value) 
  seriesNames(r) <- seriesNames(x)
  r
  }

tfSet.POSIXct <- function(value, x){ 
  if(Tobs(value) != Tobs(x)) stop("number of Tobs of observations must correspond to number of Tobs indicated by tframe.")
  r <- zoo::zoo(x, order.by = value) 
  seriesNames(r) <- seriesNames(x)
  r
  }

"seriesNames<-.zoo" <- function (x, value) 
  {if (is.matrix(x)) dimnames(x) <- list(NULL, value)
   else attr(x, "seriesNames") <- value
   x
  }

Tobs.zoo <- function(x)  NROW(x)

tfstart.zootframe <- function(x) x[1]
tfend.zootframe   <- function(x) x[length(x)]
Tobs.zootframe     <- function(x) length(x)

tfL.zoo <- function (x, p = 1) lag(x, k = -p)

tfwindow.zoo <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
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

tfExpand.zoo <- function(x, add.start = 0, add.end = 0){
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
   zoo::zoo(r, order.by = idx) 
   }

tbind.zoo <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
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


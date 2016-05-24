###############################################

#   .ts methods (for the object) and .tstframe methods (for the tframe)

###############################################
is.tframed.ts <- function(x) {TRUE}

"seriesNames<-.ts" <- function (x, value) 
  {if (is.matrix(x)) dimnames(x) <- list(NULL, value)
   else attr(x, "seriesNames") <- value
   x
  }

tframe.ts <- function(x){classed(tsp(x), c("tstframe", "tframe"))} # extractor

tfSet.tstframe <- function(value, x) {
    if (is.list(value)) { return(do.call("ts", append(list(x), value)))}
    else {
        r <- try(x <- ts(x), silent=TRUE) # vector of tsp values
        if (inherits(r, "try-error")) {attr(x, "tframe") <- value}
	else tsp(x) <- value
	return(x)
	}
    }

selectSeries.ts <- function(x, series=seqN(nseries(x))) {
  names <- seriesNames(x)
  if (is.character(series)) series <- match(names,series, nomatch=0) > 0
  if(all(0==series) | is.null(series)) r <- NULL
  else if (!is.matrix(x)) r <- x  # vector case
  else {
    r <- x[, series, drop = FALSE]
    seriesNames(r) <- names[series]
    }
  r
  }

tbind.ts <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
 {# this is used like old tsmatrix should produce a column matrix from a
  #  single vector
  nm <- seriesNames(x)
  for (z in list(...)) {
    if (!is.null(z)) {
      nm <- c(nm, seriesNames(z))
      if (!is.ts(z)) z <- ts(z,start=tfstart(z),end=tfend(z),frequency=tffrequency(z))
      x <- cbind(x, z)
      }
    }
  if (!is.matrix(x)) x <- ts(matrix(x, length(x),1),
                      start=tfstart(x), end=tfend(x), frequency=tffrequency(x))
  if (!pad.start | !pad.end)
     x <- trimNA(x, startNAs= !pad.start, endNAs= !pad.end)
  seriesNames(x) <- nm
  x
 }

tfwindow.ts <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
  {# With the default warn=T warnings will be issued if no truncation takes
   #  place because start or end is outside the range of data.
   if (!warn) 
     {opts <- options(warn = -1)
      on.exit(options(opts))
     }
   y <- window(x, start=start, end=end)
   if (is.matrix(x) && !is.matrix(y) )
      y <- tframed(matrix(y, length(y), ncol(x)), tframe(y))
   seriesNames(y) <- seriesNames(x)
   y
  }


# The next methods should work for most tstframe  tframes.
# Following are a couple that are slightly different.

tfstart.tstframe <- function(x) c(floor(x[1]+getOption("ts.eps")),
    round(1 + ((x[1]+getOption("ts.eps"))%%1)*x[3]))
#  (... further arguments, currently disregarded)

tfend.tstframe   <- function(x) c(floor(x[2]+getOption("ts.eps")),
    round(1 + ((x[2]+getOption("ts.eps"))%%1)*x[3]))
#  (... further arguments, currently disregarded)


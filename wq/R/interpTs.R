interpTs <-
function(x, type = c("linear", "series.median", "series.mean", "cycle.median", 
         "cycle.mean"), gap = NULL) {

	# Validate arguments
  gap.max <- nrow(as.matrix(x)) - 2
  if (is.null(gap))
    gap <- gap.max
	if (is.na(as.numeric(gap)) || gap < 1 || gap > gap.max)
		stop("gap must be a number between 1 and the length - 2")
  type = match.arg(type)
  if (!is.ts(x) && type %in% c("cycle.median", "cycle.mean"))
    stop("x must be a time series for these types")
 	
  # Define function for replacement by cycle
  tspx <- tsp(x)
  replaceNA <- function (x, stat) {
    x <- ts(x, start = tspx[1], frequency = tspx[3])
    x1 <- window(x, start = start(x)[1], end = c(end(x)[1], 12), extend = TRUE)
    x2 <- matrix(x1, byrow = TRUE, ncol = 12)
    stats <- apply(x2, 2, stat, na.rm = TRUE)
    indx  <- (1:length(x1))[is.na(x1)]
    x3 <- replace(x1, indx, stats[cycle(x1)[indx]])
    window(x3, start = tspx[1], end = tspx[2])
  }
  
	# Define function for vectors
	f1 <- function(x, gap, type) {
		if (sum(!is.na(x)) < 2) {
			x1 <- x
		} else {
	    x1 <- switch(type,
	      linear = na.approx(x, na.rm = FALSE),
	      series.median = ifelse(is.na(x), median(x, na.rm = TRUE), x),
	      series.mean = ifelse(is.na(x), mean(x, na.rm = TRUE), x),
        cycle.median = replaceNA(x, median),
        cycle.mean = replaceNA(x, mean)
	      )
      for (i in seq_len(length(x) - gap)) {
        seq1 <- i:(i + gap)
        if (all(is.na(x[seq1])))
          x1[seq1] <- x[seq1]
      }
		}
		x1
	}

	# Do the interpolation
	if (is.matrix(x)) {
		ans <- apply(x, 2, f1, gap, type)
	} else {
  	ans <- f1(x, gap, type)
	}
  if (is.ts(x)) {
    ans <- ts(ans, start = start(x), frequency = frequency(x))
  }
	ans
}

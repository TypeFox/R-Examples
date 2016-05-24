decompTs <-
function(x, event = TRUE, type = c("mult", "add"),
         center = c("median", "mean")) {

  # Validate input
  if (!is.ts(x) || !identical(frequency(x), 12)) {
    stop("x must be a monthly 'ts' vector")
  }
  type = match.arg(type)
  center = match.arg(center)

  # Set the time window
  startyr <- start(x)[1]
  endyr <- end(x)[1]
  x <- window(x, start = c(startyr, 1), end = c(endyr, 12), extend=TRUE)

  # Choose the arithmetic typeations, depending on type
  if (type == "mult") {
    `%/-%` <- function(x, y) x / y
    `%*+%` <- function(x, y) x * y
  } else {
    `%/-%` <- function(x, y) x - y
    `%*+%` <- function(x, y) x + y
  }

  # Choose the centering method, depending on center
  if (center == "median") {
    center <- function(x, na.rm=FALSE) median(x, na.rm=na.rm)
  } else {
    center <- function(x, na.rm=FALSE) mean(x, na.rm=na.rm)
  }

  # Long-term center
  grand <- center(x, na.rm=TRUE)

  # Annual component
  x1 <- x %/-% grand
  annual0 <- aggregate(x1, 1, center, na.rm=TRUE)
  annual1 <- as.vector(t(matrix(rep(annual0, 12), ncol=12)))
  annual <- ts(annual1, start=startyr, frequency=12)

  # Remaining components
  x2 <- x1 %/-% annual
  if (event) {
  	# Seasonal component
    seasonal0 <- matrix(x2, nrow=12)
    seasonal1 <- apply(seasonal0, 1, center, na.rm=TRUE)
    seasonal <- ts(rep(seasonal1, endyr - startyr + 1), start=startyr,
                   frequency=12)
	  # Events component
    x3 <- x2 %/-% seasonal
    # result
    ts.union(original=x, annual, seasonal, events=x3)
  } else {
    ts.union(original=x, annual, seasonal=x2)
  }
}

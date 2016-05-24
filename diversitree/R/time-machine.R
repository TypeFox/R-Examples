## TODO:
##   * accept functions, and roll back to compiled versions where
##     available.
##   * think about split functions (global variable screws with this,
##     but shifting to reassigning the pointer at initmod_*_t would
##     get around this.
make.time.machine <- function(functions, t.range, nonnegative=TRUE,
                              truncate=FALSE, k=0,
                              spline.data=NULL) {
  ## TODO: Relax this later...
  if ( !is.character(functions) )
    stop("Functions must be a character vector")
  n <- length(functions)
  nonnegative <- check.par.length(nonnegative, n)
  truncate    <- check.par.length(truncate,    n)

  ## This is pretty ugly, but should do for now.  Once we trim the old
  ## time machine code, this can be cleaned up.
  if ( "spline.t" %in% functions ) {
    if ( is.null(spline.data) )
      stop("spline data must be provided if spline.t used as a function")
    spline.data <- check.spline.data(list(t.range=t.range),
                                     spline.data)
  } else {
    if ( !is.null(spline.data) )
      warning("Ignoring spline.data -- no spline function specified")
    spline.data <- list(t=numeric(0), y=numeric(0))
  }

  tm <- new(TimeMachine, names(functions), functions, nonnegative,
            truncate, k, spline.data$t, spline.data$y)
  attr(tm, "argnames") <- time.machine.argnames(functions)
  tm
}

time.machine.argnames <- function(functions) {
  if ( !is.character(functions) )
    stop("'functions' must be a character vector")

  if ( is.null(names(functions)) ) # already true?
    stop("functions must be named")

  ## Build up argument names:
  ## Ideally we'd actually get these back from the C side too...
  info.t <- list(constant.t="c",
                 linear.t=c("c", "m"),
                 stepf.t=c("y0", "y1", "tc"),
                 sigmoid.t=c("y0", "y1", "tmid", "r"),
                 spline.t=c("y0", "y1"),
                 exp.t=c("l", "a")) ##New part added by Gustavo Burin
  if ( !all(functions %in% names(info.t)) )
    stop("Unknown functions")

  argnames <- mapply(paste, names(functions),
                     unname(info.t[functions]),
                     sep=".", SIMPLIFY=FALSE)
  is.constant <- functions == "constant.t"
  argnames[is.constant] <- names(functions)[is.constant]
  argnames <- unlist(argnames)
  names(argnames) <- NULL
  if ( any(duplicated(argnames)) )
    stop("Duplicate argument names: consider different prefixes?")
  argnames
}

check.spline.data <- function(obj, spline.data) {
  if ( !all(c("t", "y") %in% names(spline.data)) )
    stop("spline.data must contain 't' and 'y'")
  t <- spline.data$t
  y <- spline.data$y
  if ( length(t) < 3 )
    stop("Must have at least three points") # ? or 2?
  if ( length(t) != length(y) )
    stop("Lengths of t and y in spline.data must be equal")
  if ( any(is.na(t)) || any(is.na(y)) )
    stop("Neither t nor y may contain NA values")
  if ( any(duplicated(t)) )
    stop("t cannot contain duplicated values")
  t.range <- obj$t.range
  if ( min(t) > t.range[1] || max(t) < t.range[2] )
    stop(sprintf("Spline data must span time range: [%s, %s]",
                 min(t.range), max(t.range)))

  deriv <- if (is.null(spline.data$deriv))
    0L else check.integer(spline.data$deriv)
  if ( length(deriv) != 1 )
    stop("spline.data$deriv must be an scalar")

  ## Renormalise y data onto [0,1], after projecting what the actual
  ## minimum and maximum values are.
  ##
  ## TODO: This could be prone to error, especially as I've just
  ## hardcoded a "large" number of x spacing in here.  Ideally we
  ## would be able to provide this with the spline.data?  But the
  ## spline range may slightly exceed the values that are predicted
  ## analytically.
  tt <- seq(t.range[1], t.range[2], length.out=10001)
  r <- range(spline(t, y, xout=tt)$y)
  y <- (y - r[1]) / (r[2] - r[1])

  list(t=as.numeric(t), y=as.numeric(y), deriv=deriv)
}

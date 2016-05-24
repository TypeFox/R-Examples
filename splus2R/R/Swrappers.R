###################################################
# SPLUS-R Wrapper Functions for R-Fractal package
###################################################

###################################################
# UTILITY Function Definitions
#
# anyMissing
# as.rectangular
# colIds
# colMaxs
# colMedians
# colMins
# colRanges
# colStdevs
# colVars
# deparseText
# ifelse1
# is.missing
# isNumericAtomicVector
# is.rectangular
# lowerCase
# numCols
# numRows
# oldUnclass
# positions
# rmvnorm
# rowIds
# rowMaxs
# rowMins
# rowRanges
# rowStdevs
# rowVars
# stdev
# subscript2d
# upperCase
# vecnorm
# which.na
#
###################################################

###
# anyMissing
###

"anyMissing" <- function(x)
  any(is.na(unlist(x)))

###
# as.rectangular
###

"as.rectangular" <- function(x)
{
  if(is.rectangular(x))
    x
  else
    as.data.frame(x)
}

###
# colIds
# colMaxs
# colMedians
# colMins
# colRanges
# colStdevs
# colVars
# rowIds
# rowMaxs
# rowMins
# rowRanges
# rowStdevs
# rowVars
###

"colIds" <- colnames

"colMaxs" <- function(x, na.rm = FALSE, dims = 1, n = NULL)
{
  if(!identical(dims,1))
    stop("Only dims=1 is supported")
  if(!is.null(n))
    stop("Argument n is not supported")
  if(isNumericAtomicVector(x))
    return(max(x, na.rm=na.rm))
  if(!is(x,"matrix") && !is(x,"data.frame"))
    stop("Input must be a matrix or a data.frame")

  apply(x, MARGIN=2, FUN=max, na.rm=na.rm)
}

"colMedians" <- function(x, na.rm=FALSE)
{
  if(isNumericAtomicVector(x))
    return(median(x, na.rm=na.rm))
  if(!is(x,"matrix") && !is(x,"data.frame"))
    stop("Input must be a matrix or a data.frame")

  unlist(apply(x, MARGIN=2, FUN=median))
}

# this overloads the R version, which would otherwise
# choke on vectors, e.g., colMeans(1:100) would fail.
# the version below will handle this.
#"colMeans" <- function(x, na.rm=FALSE, dims=NULL)
#{
#	if(isNumericAtomicVector(x))
#	  return(mean(x, na.rm=na.rm))
#
#  # try coercion if input is the wrong class
#	if(!is(x,"matrix") && !is(x,"data.frame"))
#	  x <- as.matrix(x)
#
#	unlist(apply(x, MARGIN=2, FUN=mean))
#}

"colMins" <- function(x, na.rm = FALSE, dims = 1, n = NULL)
{
  if(!identical(dims,1))
    stop("Only dims=1 is supported")
  if(!is.null(n))
    stop("Argument n is not supported")
  if(isNumericAtomicVector(x))
    return(min(x, na.rm=na.rm))
  if(!is(x,"matrix") && !is(x,"data.frame"))
    stop("Input must be a matrix or a data.frame")

  apply(x, MARGIN=2, FUN=min, na.rm=na.rm)
}

"colRanges" <- function(x, na.rm=FALSE, dims = 1, n = NULL)
{
  if(!identical(dims,1))
    stop("Only dims=1 is supported")
  if(!is.null(n))
    stop("Argument n is not supported")
  if(isNumericAtomicVector(x))
    return(range(x, na.rm=na.rm))
  if(!is(x,"matrix") && !is(x,"data.frame"))
    stop("Input must be a matrix or a data.frame")

  apply(x, MARGIN=2, FUN=range, na.rm=na.rm)
}

#"colStdevs" <- function(x, ...) sqrt(colVars(x, ...))

"colVars" <- function(x, na.rm=FALSE, dims = 1, unbiased = TRUE,
		      SumSquares = FALSE, weights = NULL,
		      freq = NULL, n = NULL)
{
  if(!identical(dims,1))
    stop("Only dims=1 is supported")
  if(!identical(unbiased,TRUE))
    stop("Only unbiased=TRUE is supported")
  if(!identical(SumSquares,FALSE))
    stop("Only unbiased=TRUE is supported")
  if(!is.null(weights))
    stop("Argument weights is not supported")
  if(!is.null(freq))
    stop("Argument freq is not supported")
  if(!is.null(n))
    stop("Argument n is not supported")
  if(isNumericAtomicVector(x))
    return(var(x, na.rm=na.rm))

  unlist(apply(x, MARGIN=2, FUN=var, na.rm=na.rm))
}

"rowIds" <- rownames

"rowMaxs" <- function(x, na.rm = FALSE, dims = 1, n = NULL)
{
  if(!identical(dims,1))
    stop("Only dims=1 is supported")
  if(!is.null(n))
    stop("Argument n is not supported")
  if(isNumericAtomicVector(x))
    return(x)
  if(!is(x,"matrix") && !is(x,"data.frame"))
    stop("Input must be a matrix or a data.frame")

  apply(x, MARGIN=1, FUN=max, na.rm=na.rm)
}

"rowMins" <- function(x, na.rm = FALSE, dims = 1, n = NULL)
{
  if(!identical(dims,1))
    stop("Only dims=1 is supported")
  if(!is.null(n))
    stop("Argument n is not supported")
  if(isNumericAtomicVector(x))
    return(x)
  if(!is(x,"matrix") && !is(x,"data.frame"))
    stop("Input must be a matrix or a data.frame")

  apply(x, MARGIN=1, FUN=min, na.rm=na.rm)
}

"rowRanges" <- function(x, na.rm = FALSE, dims = 1, n = NULL)
{
  if(!identical(dims,1))
    stop("Only dims=1 is supported")
  if(!is.null(n))
    stop("Argument n is not supported")
  if(isNumericAtomicVector(x))
    return(rbind(x,x))
  if(!is(x,"matrix") && !is(x,"data.frame"))
    stop("Input must be a matrix or a data.frame")

  apply(x, MARGIN=1, FUN=range, na.rm=na.rm)
}

#rowStdevs <- function(x, ...) sqrt(rowVars(x, ...))

#rowVars <- function(x, na.rm = FALSE, dims = 1, unbiased = TRUE, SumSquares = FALSE,
#		    weights = NULL, freq = NULL, n = NULL)
#{
#  if(!identical(dims,1))
#    stop("Only dims=1 is supported")
#  if(!identical(unbiased,TRUE))
#    stop("Only unbiased=TRUE is supported")
#  if(!identical(SumSquares,FALSE))
#    stop("Only unbiased=TRUE is supported")
#  if(!is.null(weights))
#    stop("Argument weights is not supported")
#  if(!is.null(freq))
#    stop("Argument freq is not supported")
#  if(!is.null(n))
#    stop("Argument n is not supported")
#  if(isNumericAtomicVector(x))
#    return(NA*x)
#  if(!is(x,"matrix") && !is(x,"data.frame"))
#    stop("Input must be a matrix or a data.frame")
#
#  apply(x, MARGIN=1, FUN=var, na.rm=na.rm)
#}

###
# deparseText
###

"deparseText" <- function(expr, maxchars=30){
  # deparse the argument into a single string, with at most `maxchars' characters.
  # New lines are turned into blanks, and truncated results end in `"...."'

  full <- paste(deparse(expr), collapse=" ")
  if(nchar(full) > maxchars)
    paste(substring(full, 1, maxchars-4), "....", sep="")
  else
    full
}


###
# ifelse1
###

"ifelse1" <- function(test, x, y, ...)
{
  # if(test) return x, else return y.
  # Like ifelse(), except that test is length 1, and x or y
  # is returned as is (whatever its length).
  #     ifelse1(test, x, y)
  # is equivalent to
  #     if(test){x} else {y}.
  # This is particularly useful for assignment;
  #     answer = ifelse1(test, x, y)
  # is equivalent to
  #     if(test) answer = x else answer = y
  #
  # If more than three arguments are supplied, then y should be
  # a second test;
  #     ifelse1(test1, x, test2, u, v)
  # is equivalent to
  #     if(test){x} else if(test2) {y} else {v}
  # This may be iterated; there should be an odd number of arguments.
  if(test) x else if(missing(..1))
    y
  else ifelse1(y, ...)
}


###
# is.missing
###

"is.missing" <- function(x){
  if(length(x))
    is.na(x)
  else
    TRUE
}

###
# isNumericAtomicVector
#
# NOTE: This function is NOT defined in S-PLUS, but it is used
# often in the splus2R package and so it is defined here for
# convenience. A similar function (isVectorAtomic) is defined
# in the IFULTOOLS package, but replacing all occurrences of
# isNumericAtomicVector() with isVectorAtomic() would form a mutual
# dependence between the packages, which is not desirable.
###

"isNumericAtomicVector" <- function(x)
  is.atomic(x) && ((numRows(x) == 1 || numCols(x)) == 1) && is(x,"numeric")

###
# is.rectangular
###

"is.rectangular" <- function(x){
  # Rectangular data objects include matrices, data frames, bdFrames, atomic vectors,
 (is.matrix(x) || is.data.frame(x) || ( is.vector(x) && is.atomic(x) ))
}

###
# lowerCase
# upperCase
###

"lowerCase" <- function(x)
  casefold(x, upper=FALSE)

"upperCase" <- function(x)
  casefold(x, upper=TRUE)


###
# nDotArgs
###
"nDotArgs" <- function(...)
{
  # the number of arguments corresponding to `"..."' in the current call.
  nargs()
}

###
# numCols
# numRows
###

"numCols" <- function(x)
{
  if(is.matrix(x) || is.data.frame(x))
   ncol(x)
  else if(is.atomic(x) && is.vector(x))
   1
  else
   NULL
}

"numRows" <- function(x)
{
  if(is.matrix(x) || is.data.frame(x))
   nrow(x)
  else if((is.atomic(x) && is.vector(x)) || is(x,"signalSeries"))
   length(x)
  else
   NULL
}

###
# oldUnclass
###

"oldUnclass" <- function(x)
{
  ## the S3 version of function `unclass'; it sets `oldClass' to `NULL', rather than
  ## `class'.
  oldClass(x) <- NULL
  x
}

###
# positions
###

"positions" <- function(object)
{
  # return the positions of an ordered data object
  object@positions
}

###
# rmvnorm
###

"rmvnorm" <- function(n, mean = rep(0, d), cov = diag(d), sd, rho, d = 2)
{
  if(length(n) > 1)
    stop("n must be a single number")
  # process mean
  if(!missing(mean)) {
    if(is.matrix(mean)) {
      if(nrow(mean) != n)
	stop("mean must be vector, or matrix with n rows"
	     )
      d.mean <- ncol(mean)
    }
    else d.mean <- length(mean)
  }
  # process covariance or correlation
  method <- 1
  if(!missing(cov)) {
    if(!is.matrix(cov) || diff(dim(cov)) || any(abs(cov - t(cov)) >
						sqrt(.Machine$double.eps)))
      stop("cov must be a square symmetric matrix")
    d.cov <- nrow(cov)
    if(d.cov > 1 && any(abs(cov - diag(diag(cov))) > sqrt(.Machine$
							  double.eps)))
      method <- 3
  }
  else if(!missing(rho)) {
    if(length(rho) != 1 & length(rho) != n)
      stop("rho must have length 1 or n")
    if(any(rho < -1 | rho > 1))
      stop("rho must be between -1 and +1")
    method <- 2
    if(all(rho == 0))
      method <- 1
    d.rho <- 2
  }

  # need to add this for R, otherwise
  # (!missing(sd)) will return a different
  # answer than does S-PLUS after the if()
  # statement below has been processed
  is.missing.sd <- missing(sd)

  if(!missing(sd)) {
    if(any(sd <= 0))
      stop("Negative or zero sd found")
    if(is.matrix(sd)) {
      if(nrow(sd) != n)
	stop("number of rows of sd must equal n")
      d.sd <- ncol(sd)
    }
    else {
      d.sd <- length(sd)
      sd <- rep(sd, each = n)
    }
  }
  else sd <- 1

  # check that inferences for d match
  d.guesses <- c(d = if(!missing(d)) d, cov = if(!missing(cov)) d.cov,
		 mean = if(!missing(mean)) d.mean, rho = if(missing(cov) & !
						     missing(rho)) d.rho, sd = if(!is.missing.sd) d.sd)
  if(length(d.guesses)) {
    if(min(d.guesses) != max(d.guesses))
      stop(paste("value of d is ambiguous from arguments",
		 paste(paste(names(d.guesses), ":", d.guesses),
		       collapse = " ")))
    if(missing(d))
      d <- d.guesses[1]
  }
  # generate random numbers
  if(method == 1) {
    # Independent columns
    if(d > 1) cov <- diag(cov)
    if(any(cov < 0))
      stop("Negative variance detected")
    z <- matrix(rnorm(n * d), n, d) * rep(sqrt(cov), each = n) * sd
  }
  else if(method == 2) {
    # Use rho
    z <- matrix(rnorm(n * 2), n, 2)
    z[, 2] <- rho * z[, 1] + sqrt(1 - rho^2) * z[, 2]
    if(!missing(sd))
      z <- z * sd
  }
  else {
    # Default, multivariate method
    eS <- eigen(cov, symmetric = TRUE)
    if(any(eS$values < 0))
      stop("cov is not positive definite")
    z <- matrix(rnorm(d * n), n) %*% (sqrt(eS$values) *
				      t(eS$vectors))
    if(!missing(sd))
      z <- z * sd
  }
  if(missing(mean))
    z
  else z + (if(is.matrix(mean)) mean else rep(mean, each = n))
}

###
# stdev
###

"stdev" <- function(x, ...)
  sqrt(colVars(c(x), ...))

###
# subscript2d
###

"subscript2d" <- function(x,i,j){
  UseMethod("subscript2d")
}

subscript2dMatrix <- function(x,i,j){
  if(!missing(i) && !missing(j))
    return(x[i, j, drop = FALSE])
  if(!missing(i))
    return(x[i,  , drop = FALSE])
  if(!missing(j))
    return(x[, j, drop = FALSE])
  x[,  , drop = FALSE]
}

subscript2dDataFrame <- subscript2dMatrix


subscript2d.default <- function(x,i,j){
  # Subscript function for rectangular objects
  if(length(dim(x)) == 2){
    return(subscript2dMatrix(x, i, j))
  }
  # rest is for atomic-like vectors
  if(!missing(j)) {
    if(mode(j) == "numeric") {
      j <- j[j != 0 & j != -1]
      if(!length(j))
        return(x[0])
      if(any(j != 1))
        stop("2nd subscript out of range")
      if(length(j) > 1)
        stop("2nd subscript too long")
    }
    else if((mode(j) == "logical") && (length(j) > 1))
      stop("2nd subscript too long")
    else stop("2nd vector subscript must be numeric or logical")
  }
  if(missing(i))
    return(x[])
  len <- length(x)
  if(mode(i) == "numeric") {
    if(!length(i))
      return(x[])
    if(any(i > len | i <  - len))
      stop("1st subscript out of range")
  }
  else if((mode(i) == "logical") && (length(i) > len))
    stop("1st subscript too long")
  else if(is(i, "character") &&
	  (is.null(names(x)) || anyMissing(match(i, names(x)))))
    stop("non-matching 1st subscript")
  x[i]
}


###
# vecnorm
###

"vecnorm" <- function(x, p=2)
{
  if(is.character(p)){
    if(charmatch(p, "maximum", nomatch = 0) == 1)
      p <- Inf
    else if(charmatch(p, "euclidean", nomatch = 0) == 1)
      p <- 2
    else
      stop("improper specification of p")
  }
  if(!is.numeric(x) && !is.complex(x))
    stop("mode of x must be either numeric or complex")
  if(!is.numeric(p))
    stop("improper specification of p")
  if(p < 1)
    stop("p must be greater than or equal to 1")

  x <- ifelse1(is.numeric(x), abs(x), Mod(x))

  if(p == Inf)
    return(max(x))
  if(p == 1)
    return(sum(x))

  xmax <- max(x)
  if(!xmax)
    return(xmax)
  x <- x/xmax
  xmax * sum(x^p)^(1/p)
}


###
# which.na
###

"which.na" <- function(x)
  which(is.na(x))


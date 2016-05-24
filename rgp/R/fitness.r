## fitness.R
##   - Standard GP fitness functions and tools for creating GP fitness functions
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Mean absolute error (MAE)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The MAE between \code{x} and \code{y}.
##' @export
##' @useDynLib rgp do_mae
mae <- function(x, y) .Call(do_mae, x, y)

##' R version of Mean absolute error (MAE)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The MAE between \code{x} and \code{y}.
##' @export
r_mae <- function(x, y) mean(abs(x - y))

##' Sum squared error (SSE)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The SSE between \code{x} and \code{y}.
##' @export
##' @useDynLib rgp do_sse
sse <- function(x ,y) .Call(do_sse, x, y)

##' R version of Sum squared error (SSE)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The SSE between \code{x} and \code{y}.
##' @export
r_sse <- function(x ,y) sum((x -y) ^ 2)

##' Mean squared error (MSE)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The MSE between \code{x} and \code{y}.
##' @export
mse <- function(x, y) .Call(do_sse, x, y) / length(x)

##' Root mean squared error (RMSE)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The RMSE between \code{x} and \code{y}.
##' @export
rmse <- function(x, y) sqrt(mse(x, y))

##' Scaled sum squared error (sSSE)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The sSSE between \code{x} and \code{y}.
##' @export
##' @useDynLib rgp do_ssse
ssse <- function(x, y) .Call(do_ssse, x, y)

##' R version of Scaled sum squared error (sSSE)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The sSSE between \code{x} and \code{y}.
##' @export
r_ssse <- function(x, y) {
	if (length(y) == 1) y <- rep(y, length(x))
	b = cov(x, y) / var(y)
	a = mean(x) - b * mean(y)
	sum((x - (a + b * y)) ^ 2)
}

##' Coefficient of determination (R^2)
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The coefficient of determination (R^2) between \code{x} and \code{y}.
##' @export
rsquared <- function(x, y) {
  ssTot <- sum((x - mean(x)) ^ 2)
  ssRes <- sum((x - y) ^ 2)
  1 - ssRes / ssTot
}

##' Normalize a vector into the interval [0, 1]
##'
##' @param x The vector to normalize, so that each element lies in the
##'   interval [0, 1].
##' @return The normalized vector.
##' @export
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

##' Normalized mean squared error (NMSE)
##'
##' Calculates the MSE between vectors after normalizing them into the
##' interval [0, 1].
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The NMSE between \code{x} and \code{y}.
##' @export
nmse <- function(x, y) mse(normalize(x), normalize(y))

##' Scaled mean squared error (SMSE)
##'
##' Calculates the MSE between vectors after scaling them.
##' Beware that this error measure is invariant to scaling with negative constants, i.e.
##' the multiplicative inverse of the true functions also receives an error of 0. 
##' See \url{http://www2.cs.uidaho.edu/~cs472_572/f11/scaledsymbolicRegression.pdf} for
##' details.
##'
##' @param x A numeric vector or list.
##' @param y A numeric vector or list.
##' @return The NMSE between \code{x} and \code{y}.
##' @export
smse <- function(x, y) (1 / length(x)) * ssse(x, y)

##' Symbolic squared error function (SE)
##'
##' Given two functions \code{f} and \code{g}, returns a function whose body
##' is the symbolic representation of the squared error between \code{f} and
##' \code{g}, i.e. \code{function(x) (f(x) - g(x))^2}.
##'
##' @param f An R function.
##' @param g An R function with the same formal arguments as \code{f}.
##' @return A function representing the squared error between \code{f} and \code{g}. 
##' @export
seSymbolicFunction <- function(f, g) {
  newf <- new.function()
  formals(newf) <- formals(f)
  body(newf) <- call("^", call("-", body(f), body(g)), 2)
  newf
}

##' Symbolic squared error (SE)
##'
##' Given to functions \code{f} and \code{g}, returns the area the squared
##' differences between \code{f} and \code{g} in the integration limits
##' \code{lower} and \code{upper}.
##'
##' @param f An R function.
##' @param g An R function with the same formal arguments as \code{f}.
##' @param lower The lower limit of integraion.
##' @param upper The upper limit of integraion.
##' @param subdivisions The maximum number of subintervals for numeric integration. 
##' @return The area of the squared differences between \code{f} and \code{g}, or
##'   \code{Inf} if integration is not possible in the limits given. 
##' @export
seSymbolic <- function(f, g, lower, upper, subdivisions = 100) {
  seFunction <- seSymbolicFunction(f, g)
  tryCatch(integrate(seFunction, lower = lower, upper = upper, subdivisions = subdivisions)$value,
           error = function(error) Inf)
           #error = function(error) { print(error); Inf }) # DEBUG
}

##' Create a fitness function based on symbolic squared error (SE) 
##'
##' Creates a fitness function that calculates the squared error of
##' an individual with respect to a reference function \code{func}.
##' When an \code{indsizelimit} is given, individuals exceeding this
##' limit will receive a fitness of \code{Inf}.
##'
##' @param func The reference function.
##' @param lower The lower limit of integraion.
##' @param upper The upper limit of integraion.
##' @param subdivisions The maximum number of subintervals for numeric integration. 
##' @param indsizelimit Individuals exceeding this size limit will get
##'   a fitness of \code{Inf}.
##' @return A fitness function based on the reference function \code{func}.
##' @export
makeSeSymbolicFitnessFunction <- function(func, lower, upper, subdivisions = 100,
                                          indsizelimit = NA) {
  function(ind) {
    errorind <- seSymbolic(ind, func,
                           lower = lower, upper = upper,
                           subdivisions = subdivisions)
  	if (!is.na(indsizelimit) && funcSize(ind) > indsizelimit)
  	  Inf # ind size limit exceeded
  	else errorind
  }
}

##' Create a fitness function from a reference function of one variable
##'
##' Creates a fitness function that calculates an error measure with
##' respect to an arbitrary reference function of one variable on the
##' sequence of fitness cases \code{seq(from, to, length = steps)}.
##' When an \code{indsizelimit} is given, individuals exceeding this
##' limit will receive a fitness of \code{Inf}.
##'
##' @param func The reference function.
##' @param from The start of the sequence of fitness cases.
##' @param to The end of the sequence of fitness cases.
##' @param steps The number of steps in the sequence of fitness cases.
##' @param errorMeasure A function to use as an error measure, defaults to RMSE.
##' @param indsizelimit Individuals exceeding this size limit will get
##'   a fitness of \code{Inf}.
##' @return A fitness function based on the reference function \code{func}.
##' @export
makeFunctionFitnessFunction <- function(func, from = -1, to = 1, steps = 128,
                                        errorMeasure = rmse, indsizelimit = NA) {
  xs <- seq(from, to, length = steps)
  ystarget <- func(xs)
  function(ind) {
    ysind <- ind(xs) # vectorized fitness-case evaluation
  	errorind <- errorMeasure(ystarget, ysind)
  	if (!is.na(indsizelimit) && funcSize(ind) > indsizelimit)
  	  Inf # ind size limit exceeded
  	else if (is.na(errorind) || is.nan(errorind))
  	  Inf # error value is NA or NaN
  	else errorind
  }
}

##' Create a fitness function from a n-ary reference function
##'
##' Creates a fitness function that calculates an error measure with
##' respect to an arbitrary n-ary reference function based sample points
##' generated by a given \code{designFunction}.  
##' When an \code{indsizelimit} is given, individuals exceeding this
##' limit will receive a fitness of \code{Inf}.
##'
##' @param func The reference function. Its single argument must be numeric
##'   vector of length \code{dim} and it must return a scalar numeric.
##' @param dim The dimension of the reference function.
##' @param designFunction A function to generate sample points. Its first
##'   argument must be \code{dim}. Defaults to \code{\link{gridDesign}}.
##' @param errorMeasure A function to use as an error measure, defaults to RMSE.
##' @param indsizelimit Individuals exceeding this size limit will get
##'   a fitness of \code{Inf}.
##' @param ... Additional arguments to the \code{designFunction}.
##' @return A fitness function based on the reference function \code{func}.
##'
##' @seealso \code{\link{latinHypercubeDesign}}, \code{\link{gridDesign}},
##' @export
makeNaryFunctionFitnessFunction <- function(func, dim, designFunction = gridDesign,
                                            errorMeasure = rmse, indsizelimit = NA, ...) {
  xs <- designFunction(dim, ...) 
  ystarget <- apply(xs, 1, func) # apply func row-wise to xs
  function(ind) {
    ysind <- apply(xs, 1, ind) # row-wise fitness-case evaluation
  	errorind <- errorMeasure(ystarget, ysind)
  	if (!is.na(indsizelimit) && funcSize(ind) > indsizelimit)
  	  Inf # ind size limit exceeded
  	else if (is.na(errorind) || is.nan(errorind))
  	  Inf # error value is NA or NaN
  	else errorind
  }
}

##' Tools for manipulating boolean functions
##'
##' \code{integerToBoolean} converts a scalar positive integer (or zero) to its
##' binary representation as list of logicals.
##' \code{booleanFunctionVector} returns the boolean vector of result values of
##' \code{f}, given a boolean function \code{f}.
##' \code{numberOfDifferentBits} given two lists of booleans of equal length,
##' returns the number of differing bits.
##' \code{makeBooleanFitnessFunction} given a boolean target function, returns
##' a fitness function that returns the number of different places in the output
##' of a given boolean function and the target function.
##'
##' @param i A scalar positive integer. 
##' @param width The with of the logical vector to return.
##' @param f A boolean function.
##' @param a A list of booleans.
##' @param b A list of booleans.
##' @param targetFunction A boolean function.
##' @return The function result as described above.
##'
##' @rdname booleanFunctionManipulation 
##' @export
integerToLogicals <- function(i, width = floor(log(base = 2, i) + 1)) {
  k <- i
  result <- list()
  l <- 0
  while (k != 0) {
    bit <- as.logical(k %% 2)
    result <- c(bit, result)
    k <- k %/% 2
    l <- l + 1
  }
  if (width - l < 0) stop("integerToLogicals: width to small")
  c(replicate(width - l, FALSE), result) # prefix with zeros
}

##' @rdname booleanFunctionManipulation 
##' @export
booleanFunctionAsList <- function(f) {
  fArity <- length(formals(f))
  inputs <- Map(function(i) integerToLogicals(i, width = fArity), 0:(2 ^ fArity - 1))
  Map(function(input) do.call(f, input), inputs)
}

##' @rdname booleanFunctionManipulation 
##' @export
numberOfDifferentBits <- function(a, b) {
  Reduce(`+`, Map(function(ai, bi) if (ai != bi) 1 else 0, a, b))
}

##' @rdname booleanFunctionManipulation 
##' @export
makeBooleanFitnessFunction <- function(targetFunction) {
  targetFunctionList <- booleanFunctionAsList(targetFunction)
  function(f) {
    fList <- booleanFunctionAsList(f)
    numberOfDifferentBits(fList, targetFunctionList)
  }
}


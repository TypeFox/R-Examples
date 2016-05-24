
# the prototype of the function defined is required in this order 
# for consistency across methods
# for the default method the prototype 
# (x, add = TRUE, morder, ...)
# would be more appropriate since "morder" is not expected to be used, 
# but for ".Arima" (used by the functions of the package) the prototype 
# (x, morder, add = TRUE, ...)
# is better since the argument for which a default value can be assigned 
# comes after the other arguments

coefs2poly <- function(x, morder, add = TRUE, ...) UseMethod("coefs2poly")

coefs2poly.Arima <- function(x, morder = NULL, add = TRUE, ...)
{
  coefs2poly(coef(x), x$arma, add = add, ...)
}

coefs2poly.default <- function(x, morder, add = TRUE, ...)
{
  if (missing(morder))
    stop("model order is not specified")

  arcoefs <- x[grep("^ar\\d+$", names(x))]
  #stopifnot(length(arcoefs) == morder[1])

  ind <- grep("^sar\\d+$", names(x))
  #stopifnot(length(ind) == morder[3])
  lenind <- length(ind)
  if (lenind > 0)
  {
    slenind <- morder[5] * lenind
    sarcoefs <- rep(0, slenind)
    sarcoefs[seq.int(morder[5], slenind, morder[5])] <- x[ind]
    if (length(arcoefs) > 0) {
      arcoefs <- -coef(polynomial(c(1, -arcoefs)) * polynomial(c(1, -sarcoefs)))[-1]
    } else
      arcoefs <- sarcoefs
  }

  macoefs <- x[grep("^ma\\d+$", names(x))]
  #stopifnot(length(macoefs) == morder[2])

  ind <- grep("^sma\\d+$", names(x))
  #stopifnot(length(ind) == morder[4])
  lenind <- length(ind)
  if (lenind > 0)
  {
    slenind <- morder[5] * lenind
    smacoefs <- rep(0, slenind)
    smacoefs[seq.int(morder[5], slenind, morder[5])] <- x[ind]
    if (length(macoefs) > 0) {
      macoefs <- coef(polynomial(c(1, macoefs)) * polynomial(c(1, smacoefs)))[-1]
    } else
      macoefs <- smacoefs
  }

  if (isTRUE(add) && morder[6] == 1)
  {
    #multiply arcoefs by (1 - L)
    if (length(arcoefs) == 0) {
      # do not do minus 1 since it is moved to the right hand side of the equation
      arcoefs <- 1
    } else {
      # multiply (1 - arcoefs(L)) by (1 - L)
      # based on package "polynom"
      # it is assumed that 'arcoefs' are related to consecutive lags
      tmp <- outer(c(1, -arcoefs), c(1, -1))
      arcoefs <- as.vector(tapply(tmp, row(tmp) + col(tmp), sum))
      # remove L^0 = 1 and move to the right hand side of the equation
      arcoefs <- -arcoefs[-1]
    }
  }

  if (isTRUE(add) && morder[7] == 1)
  {
    # multiply arcoefs by (1 - L)
    if (length(arcoefs) == 0) {
      arcoefs <- c(rep(0, morder[5] - 1), 1)
    } else {
      # multiply (1 - arcoefs(L)) by (1 - L)
      # based on package "polynom"
      # it is assumed that 'arcoefs' are related to consecutive lags
      tmp <- outer(c(1, -arcoefs), c(1, rep(0, morder[5] - 1), -1))
      arcoefs <- as.vector(tapply(tmp, row(tmp) + col(tmp), sum))
      # remove L^0 and move to the right hand side of the equation
      arcoefs <- -arcoefs[-1]
    }
  }

##FIXME TODO
  if (morder[6] > 1)
    stop("unsupported number of regular differences.")
  if (morder[7] > 1)
    stop("unsupported number of seasonal differences.")

  structure(list(arcoefs = arcoefs, macoefs = macoefs), class = "ArimaPars")
}

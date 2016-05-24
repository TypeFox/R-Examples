is.stepMax <- function(x, tol = .Machine$double.eps^0.5){  # test for stepMax
  if(!is.numeric(x)){stop("The stepMax you entered is not numeric!")}
  if(!is.na(x[2])){stop("The stepMax you entered is a vector!")}
  ret <- (abs(x - trunc(x)) < tol)
  if(!ret){stop("Your stepMax is not a integer!")}
  ret <- (x > 0)
  if (!ret){stop("The stepMax is negative or zero!")}
}
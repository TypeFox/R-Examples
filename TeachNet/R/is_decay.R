is.decay <- function(x, tol = .Machine$double.eps^0.5){  # test for decay parameter
  if(!is.numeric(x)){stop("The decay you entered is not numeric!")}
  if(!is.na(x[2])){stop("The decay you entered is a vector!")}
  if(!(x>=0)){stop("The decay you entered is less than zero!")}
}
is.sample <- function(x){ # test for sameSample parameter
  if(!is.logical(x)){stop("The sameSample you entered is not a logical value!")}
  if(!is.na(x[2])){stop("The sameSample you entered is a vector!")}
}
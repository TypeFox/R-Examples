is.numberOfNeurons <- function(x, tol = .Machine$double.eps^0.5){  # test for hidden.structure
  if(!is.numeric(x)){stop("The hidden.structure you entered is not a numeric vector!")}
  if(!is.na(x[3])){stop("The hidden.structure you entered is not a valid structure!")}
  if(!((x[1]==-1) && is.na(x[2]))){
    ret <- (abs(x - trunc(x)) < tol)
    if (!is.na(unique(ret)[2]) || !unique((ret[1]))){stop("There is a non integer number in hidden.structure!")}
    ret <- (x >= 0)
    if (!is.na(unique(ret)[2]) || !unique((ret[1]))){stop("There is a unallowed negative value in hidden.structure!")}
  }
}
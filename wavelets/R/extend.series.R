extend.series <- function(X, method="reflection", length="double", n, j){

  # error checking
  if(is.na(match(class(X), c("numeric", "ts", "matrix", "data.frame"))))
     stop("Invalid argument: 'X' must be of class 'numeric','ts', 'matrix', or 'data.frame'.")
  if(is.na(match(method, c("periodic","reflection","zeros","mean","reflection.inverse"))))
    stop("Invalid argument value for 'method'")
  if(is.na(match(length, c("arbitrary","powerof2","double"))))
    stop("Invalid argument value for 'length'")
  if(length == "arbitrary"){
    if(missing(n)) stop("Unspecified argument: argument 'n' must be specified when length='arbitrary'.")
    else if(round(n) != n) stop("Invalid argument: 'n' must be an integer value.")
  }
  if(length == "powerof2"){
    if(missing(j)) stop("Unspecified argument: argument 'j' must be specified when length='powerof2'.")
    if(round(j) != j) stop("Invalid argument: 'j' must be an integer value.")
  }
  
  # store the old class for output
  class.X <- class(X)
  if(class.X != "matrix"){
    attr.X <- attributes(X)
    X <- as.matrix(X)
  }
 dim.X <- dim(X)
 if((dim.X[1] == 1) & (dim.X[2] > 1)){
    X <- t(X)
    N <- dim.X[2]
  } else N <- dim.X[1]

  # determine final length 'n' of series after extension
  if(length == "arbitrary"){
    if(n <= N) stop("Invalid argument: 'n' must be greater than length of series when length='arbitrary'.")
  } else if(length == "powerof2"){
    k <- N/(2^j)
    if(round(k) == k) stop("Invalid argument: length of series should not be multiple of 2^j when length='powerof2'.")
    else n <- ceiling(k)*2^j
  } else if(length == "double") n <- 2*N

  # extend the series to length 'n'
  if(method == "periodic") X <- apply(X, 2, rep, length=n)
  if(method == "reflection") X <- apply(X, 2, function(x,n){rep(c(x,rev(x)),length=n)}, n=n)
  if(method == "zeros") X <- apply(X, 2, function(x,n){c(x,rep(0,length=n-N))}, n=n)
  if(method == "mean") X <- apply(X, 2, function(x,n){c(x,rep(mean(x),length=n-N))}, n=n)
  if(method == "reflection.inverse") X <- apply(X, 2, function(x,n,N){rep(c(x,2*x[N]-rev(x)),length=n)}, n=n, N=N)

  # return X in orginal form
  if(class.X != "matrix"){
    if(class.X == "data.frame") X <- as.data.frame(X)
    attributes(X) <- attr.X
  }  
  return(X)
}

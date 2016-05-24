modwt <- function(X, filter="la8", n.levels, boundary="periodic", fast=TRUE){

  # error checking
  if(is.na(match(class(X)[1], c("numeric", "ts", "mts", "matrix", "data.frame"))))
    stop("Invalid argument: 'X' must be of class 'numeric','ts', 'mts', 'matrix', or 'data.frame'.")
  if(is.na(match(class(filter), c("character", "wt.filter", "numeric", "integer"))))
    stop("Invalid argument: 'filter' must be of class 'character','wt.filter','numeric', or 'integer'.")
  if(class(filter) == "wt.filter"){
    if(filter@transform == "dwt")
      stop("Invalid argument: 'transform' element of agrument 'filter' must be equal to 'modwt'.")
  }
  if(is.na(match(boundary, c("periodic","reflection"))))
    stop("Invalid argument value for 'boundary'")
  if(!is.null(dim(X))){
    if(dim(X)[1] == 1){
      stop("Unacceptable dimensions: number of observations for columns of X must be greater than 1.")
    }
  }

  # get wavelet coeficients and length
  if(class(filter) == "character" | class(filter) == "numeric" | class(filter) == "integer")
    filter <- wt.filter(filter, modwt=TRUE)
  L <- filter@L

  # convert X to a matrix
  class.X <- class(X)[1]
  attr.X <- attributes(X)
  if(is.null(attr.X)) attr.X <- list()
  if(class.X == "mts") attributes(X)[3:4] <- NULL else X <- as.matrix(X)
  dim.X <- dim(X)
  N <- dim(X)[1]
  n.series <- dim(X)[2]

  # determine the level of decomposition
  if(missing(n.levels)){
    J <- as.integer(floor(log(((N-1)/(L-1))+1)/log(2)))
  } else if(!is.numeric(n.levels) | (round(n.levels) != n.levels)){
    stop("Invalid argument value: 'n.levels' must be an integer value")
  } else J <- as.integer(n.levels)

  # reflect X for reflection method
  if(boundary == "reflection"){
    X <- extend.series(X)
    N <- 2*N
  }

  # initialize variables for pyramid algorithm
  Vj <- X
  n.boundary <- rep(NA, J)
  W.coefs <- as.list(rep(NA, length=J))
  names(W.coefs) <- lapply(1:J, function(j, x){names(x)[j] <- paste("W",j,sep="")}, x=W.coefs)
  V.coefs <- as.list(rep(NA, length=J))
  names(V.coefs) <- lapply(1:J, function(j, x){names(x)[j] <- paste("V",j,sep="")}, x=V.coefs)

  # implement the pyramid algorithm
  for(j in 1:J){
    if(fast){
      Wout <- Vout <- rep(0, length=dim(Vj)[1])
      analysis <- sapply(1:n.series,
                         function(i,v,f,j,Wout,Vout){
                           out <- .C("modwt_forward", as.double(v[,i]),
                                     as.integer(length(v[,i])),
                                     as.integer(j), as.double(f@h),
                                     as.double(f@g), as.integer(f@L),
                                     as.double(Wout), as.double(Vout),
                                     PACKAGE="wavelets")
                           return(c(out[[7]],out[[8]]))
                         }, v=Vj, f=filter, j=j, Wout=Wout, Vout=Vout)
    } else {
      analysis <- sapply(1:n.series,
                         function(i,v,f,j){
                           out <- modwt.forward(v[,i],f,j)
                           return(c(out$W,out$V))
                         }, v=Vj, f=filter, j=j)
    }
    Wj <- matrix(analysis[1:N,], ncol=n.series)
    Vj <- matrix(analysis[(N+1):(2*N),], ncol=n.series)
    W.coefs[[j]] <- Wj
    V.coefs[[j]] <- Vj
    Lj <- (2^j-1)*(L-1)+1
    n.boundary[j] <- min(Lj,N)
  }

  # creat modwt object for ouput
  modwt <- new("modwt",
               W = W.coefs,
               V = V.coefs,
               filter = filter,
               level = J,
               n.boundary = n.boundary,
               boundary = boundary,
               series = X,
               class.X = class.X,
               attr.X = attr.X,
               aligned = FALSE,
               coe = FALSE)

  return(modwt)
}

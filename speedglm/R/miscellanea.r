
control<-function(B, symmetric = TRUE, tol.values = 1e-07, tol.vectors = 1e-07,
                  out.B = TRUE, method = c("eigen","Cholesky")){
  method <- match.arg(method)
  if (!(method %in% c("eigen","Cholesky"))) 
    stop("method not valid or not implemented")
  if (method=="eigen"){
    n <- ncol(B)
    sa <- 1:n
    nok <- NULL
    auto <- eigen(B, symmetric, only.values = TRUE)
    totcoll <- sum(abs(auto$values) < tol.values)
    ncoll <- totcoll
    rank <- n - ncoll
    i <- 1
    while (ncoll != 0) {
      auto <- eigen(B, symmetric)
      j <- as.matrix(abs(auto$vectors[, n]) < tol.vectors)
      coll <- which(!j)
      coll <- coll[length(coll)]
      B <- B[-coll, -coll]
      nok[i] <- coll
      ncoll <- sum(abs(auto$values) < tol.values) - 1
      n <- ncol(B)
      i <- i + 1
    }
    ok <- if (!is.null(nok))
      sa[-nok] else sa
  }
  if (method=="Cholesky"){
    A <- chol(B, pivot = TRUE)
    pivot <- attributes(A)$"pivot"
    rank <- attributes(A)$"rank"
    ok <- sort(pivot[1:rank])
    nok <- if (rank<length(pivot)) pivot[(rank+1):length(pivot)]  else NULL
    B <- B[ok,ok]
  }
  rval <- if (out.B) list(XTX = B, rank = rank, pivot = c(ok, nok)) else
    list(rank = rank, pivot =  c(ok, nok))
  
  rval
}


cp <- function(X, w = NULL, row.chunk = NULL, sparse = FALSE){
  if (sparse) {
    if (!is.null(w)) X<-sqrt(w) * X
    if (class(X) != "dgCMatrix")
      X <- as(X, "dgCMatrix")
    new.B <- crossprod(X)
  }
  else {
    if (is.null(row.chunk)) {
      new.B <- if (is.null(w))
        crossprod(X) else crossprod(sqrt(w) * X)
    } else {
      if (row.chunk >= (nrow(X) - 1))
        row.chunk <- nrow(X)
      if (row.chunk <= 1)
        row.chunk <- 2
      mod <- nrow(X)%%row.chunk
      last.block <- (mod > 0)
      G <- nrow(X)%/%row.chunk - (last.block * (2 * row.chunk <=
                                                  nrow(X)))
      new.B <- matrix(0, ncol(X), ncol(X))
      a <- row.chunk
      j <- 1
      if (is.null(w)) {
        for (i in 1:G) {
          B <- crossprod(X[j:(i * a), ])
          new.B <- new.B + B
          j <- j + a
        }
        if (mod > 0)
          new.B <- new.B + crossprod(X[j:nrow(X), ])
      } else {
        for (i in 1:G) {
          B <- crossprod(sqrt(w[j:(i * a)]) * X[j:(i *
                                                     a), ])
          new.B <- new.B + B
          j <- j + a
        }
        if (mod > 0)
          new.B <- new.B + crossprod(sqrt(w[j:nrow(X)]) *
                                       X[j:nrow(X), ])
      }
    }
  }
  as(new.B, "matrix")
}

is.sparse <- function(X,sparselim=.9,camp=.05){
  if (prod(dim(X))<500) camp <- 1
  subX<-sample(X,round((nrow(X)*ncol(X)*camp),digits=0),replace=FALSE)
  p<-sum(subX==0)/prod(length(subX))
  if (p > sparselim) sparse <- TRUE else sparse <- FALSE
  attr(sparse,"proportion of zeros in the sample")<-p
  sparse
}


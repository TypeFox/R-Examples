uniquePrim <- function(x){    #OK
  unique.default(x)
}

unlistPrim <- function(l){    #OK
  unlist(l, use.names=FALSE)
}


setdiffPrim <- function (x, y){
  unique.default(if (length(x) || length(y))
                 x[match(x, y, 0L) == 0L]
  else x)
}

intersectPrim <- function (x, y){
  unique.default(y[match(x, y, 0L)])
}


outerPrim <- function(X,Y){
  nX  <- length(X)
  nY  <- length(Y)
  Y   <- rep(Y, rep.int(length(X), length(Y)))
  X   <- rep(X, times = ceiling(length(Y)/length(X)))
  ans <-X*Y
  dim(ans)<-c(nX,nY)
  ans
}


matchPrim<-function(x,table){ # Never used
  match(x,table)
}

aveslope <- function(obj,a=0,b=1){
  if(a>=b) stop("Argument 'b' must be greater than argument 'a'. Sorry.")
  a <- (a-obj$xmin)/obj$xrng
  b <- (b-obj$xmin)/obj$xrng
  u <- as.matrix(obj$postdraws[,grep("u.[0-9]+",names(obj$postdraws))])
  M <- obj$m
  k <- 0:(M-1)
  vec <- exp(lchoose(M-1,k) + lbeta(k+1,M-k) + log(pbeta(b,k+1,M-k)-pbeta(a,k+1,M-k)))
  S <- u%*%vec
  return(S)
  }

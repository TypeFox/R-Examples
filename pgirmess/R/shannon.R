shannon <- function (vect,base=2){
  vect <- as.numeric (vect)
  if (all (vect <= 0)) return (c(H=NA, J=NA))
  vect <- vect/sum(vect)
    vect <- vect * log(vect, base)
    h <- sum(vect[is.finite(vect)])
    hmax <- log(1/length(vect), base)
    res<-c(H = -h, J = h/hmax)
    attributes(res)$baseLog<-base
    res
  }
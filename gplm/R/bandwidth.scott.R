"bandwidth.scott" <- function(x,kernel="biweight",product=TRUE){

  x <- as.matrix(x)  ## nxd
  d <- ncol(x)
  n <- nrow(x)

  if (kernel=="gaussian"||kernel=="normal"){  ## gaussian kernel
    fac <- 1
  }else{  ## non-gaussian kernel
    fac <- kernel.constants(kernel=kernel,d=d)$d0 / (1/(2*sqrt(pi))^d) ^(1/(d+4))
    ## computes the Marron & Nolan type rescaling factor
  }
  s <- rbind( apply(x,2,sd), diff( apply(x,2,quantile,c(0.25,0.75)) )/1.349 )
  s <- apply(s,2,min)
  bandwidth <- fac *s* n^(-1/(d+4)) ## Scott's ROT (see Scott, p. 152)

  return(bandwidth)
}

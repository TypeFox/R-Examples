"kernel.function" <- function(u,kernel="biweight",product=TRUE){

  if (kernel=="triangular"){ kernel <- "triangle" }
  if (kernel=="rectangle" || kernel=="rectangular"){ kernel <- "uniform" }
  if (kernel=="quartic"){ kernel <- "biweight" }
  if (kernel=="normal"){  kernel <- "gaussian" }

  kernel.names <- c("triangle","uniform","epanechnikov","biweight",
                    "triweight","gaussian")

  c1 <- c(1,0.5,0.75,0.9375,1.09375,NA) ## cf. Wand & Jones, p175
  pp <- c(1,2,2,2,2,0)
  qq <- c(1,0,1,2,3,NA)
  names(c1) <- names(pp) <- names(qq) <- kernel.names

  if (is.null(dim(u))){
    d <- 1
    u <- matrix(u,length(u),1)
  }else{
    u <- as.matrix(u)
    d <- ncol(u)
  }
  p <- pp[kernel]
  q <- qq[kernel]
  ##print(paste("d,p,q:", d,p,q))

  volume.d <- pi^(d/2)/gamma(d/2+1)  ## volume of d-dim. unit sphere
  r1 <- c(d+1,1,(d+2),(d+2)*(d+4),(d+2)*(d+4)*(d+6),NA)
  r2 <- c(1,1,2,8,48,NA)
  names(r1) <- names(r2) <- kernel.names

  if (p >0){
    if (product){
      x <- 1-sqrt(u*u)^p
      c <- c1[kernel] 
      k <- (c^d) * apply(x^q,1,prod) * apply(x>=0,1,prod)
    }else{
      x <- 1-sqrt(rowSums(u*u))^p
      c <- r1[kernel] / (r2[kernel]*volume.d)
      k <- c * x^q * (x>=0)
    }
  ##print(paste("c:",c))
  }else{
    k <- apply(dnorm(u),1,prod)
  }
  return(k)
}

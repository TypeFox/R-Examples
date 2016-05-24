"kernel.constants" <- function(kernel="biweight",d=1,product=TRUE){

  if (kernel=="triangular"){ kernel <- "triangle" }
  if (kernel=="rectangle" || kernel=="rectangular"){ kernel <- "uniform" }
  if (kernel=="quartic"){ kernel <- "biweight" }
  if (kernel=="normal"){  kernel <- "gaussian" }

  kernel.names <- c("triangle","uniform","epanechnikov","biweight",
                    "triweight","gaussian")

  m2 <- n2 <- d0 <- NA
  volume.d <- pi^(d/2)/gamma(d/2+1)  ## volume of d-dim. unit sphere

  mm2 <- c(1/6, 1/3, 0.2, 1/7, 1/9, 1)
  nn2 <- c(2/3, 0.5, 0.6, 5/7, 350/429, 1/(2*sqrt(pi)))
                              ## cf. Wand & Jones, p. 176
  pp <- c(1,2,2,2,2,0)
  qq <- c(1,0,1,2,3,NA)
  names(mm2) <- names(nn2) <- names(pp) <- names(qq) <- kernel.names

  p <- pp[kernel]
  q <- qq[kernel]

  if (product || p==0){
    m2 <-  mm2[kernel]
    n2 <-  nn2[kernel]
    d0 <- (n2^d/(m2^2))^(1/(d+4))
  }else{
    if (p==1){
      m2 <- (d+1)/((d+2)*(d+3))
      n2 <- 2*(d+1)/((d+2)*volume.d)
    }
    if (p==2){
      m2 <- 1/(2*q+d+2)
      if (q==0){ n2 <- 1/volume.d }  ## uniform
      if (q==1){ n2 <- (4+2*d)/((4+d)*volume.d)}  ## epanechnikov
      if (q==2){ n2 <- c(5/7, 9/(5*pi), 35/(22*pi), 24/(5*pi^2), 2835/(572*pi^2),120/(7*pi^3))[d]}  ## biweight
      if (q==3){ n2 <- c(350/429, 16/(7*pi), 315/(143*pi), 50/(7*pi^2), 3465/(442*pi^2),200/(7*pi^3))[d]}  ## triweight
    }
    d0 <- (n2/(m2^2))^(1/(d+4))
  }

  return(list(m2=m2,n2=n2,d0=d0))
}

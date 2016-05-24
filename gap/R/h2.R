k <- function(r,N,adjust=TRUE)
{
  r2 <- r^2
  n <- N-1
  k1 <- ifelse(adjust,r-r*(1-r2)/2/n,r)
  k2 <- (1-r2)^2/n*(1+11*r2/2/n)
  invisible(c(k1,k2))
}

h2 <- function(mzDat=NULL,dzDat=NULL,rmz=NULL,rdz=NULL,nmz=NULL,ndz=NULL,selV=NULL)
{
  if(!is.null(mzDat))
  {
    r1 <- cor(mzDat[selV[1]],mzDat[selV[2]], use="complete")
    n1 <- length(!is.na(c(mzDat[selV[1]],mzDat[selV[2]])))
  } else {
    if(is.null(rmz)|is.null(nmz)) stop("Either raw data or correlation/sample size is neeeded")
    r1 <- rmz
    n1 <- nmz
  }
  if(!is.null(dzDat))
  {
    r2 <- cor(dzDat[selV[1]],dzDat[selV[2]], use="complete")
    n2 <- length(!is.na(c(dzDat[selV[1]],dzDat[selV[2]])))
  } else {
    if(is.null(rdz)|is.null(ndz)) stop("Either raw data or correlation/sample size is neeeded")
    r2 <- rdz
    n2 <- ndz
  }
  kmz <- k(r1,n1)
  k1mz <- kmz[1]
  k2mz <- kmz[2]
  kdz <- k(r2,n2)
  k1dz <- kdz[1]
  k2dz <- kdz[2]
  h2 <- 2 * (k1mz - k1dz)
  vh <- 4 * (k2mz + k2dz)
  c2 <- 2 * k1dz - k1mz
  vc <- 4 * k2dz + k2mz
  e2 <- 1 - k1mz
  ve <- k2mz
  ACEr_est <- as.matrix(c(h2,c2,e2,vh,vc,ve))
  rownames(ACEr_est) <- c("h2","c2","e2","vh","vc","ve")
  invisible(t(ACEr_est))
}

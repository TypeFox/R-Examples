makew <- function(shpfile=NULL,coormat=NULL,method="queen",
  knum=10,ringdist=.25,kern="tcub",window=.10,eigenvalues=FALSE) {
  library(spdep)

  if (!identical(shpfile,NULL)) {coormat <- coordinates(shpfile) }
  n = nrow(coormat)
  wmat <- array(0,dim=c(n,n))
  eigvar <- NULL

  contig = TRUE
  if (method=="rook") {contig = FALSE}

  if (method=="queen"|method=="rook") {
    neighbors <- poly2nb(shpfile,queen=contig)
    wmat <- nb2mat(neighbors,zero.policy=TRUE,style="B")
  }

  if (method=="knear") {
    kmat <- knearneigh(coormat,k=knum,longlat=TRUE)
    neighbors <- knn2nb(kmat)
    wmat <- nb2mat(neighbors,zero.policy=TRUE,style="B")
  }

  if (method=="ring") {
    for (i in seq(2,n)) {
      ni = i-1
      dist <- geodistance(longvar=coormat[1:ni,1],latvar=coormat[1:ni,2],lotarget=coormat[i,1],latarget=coormat[i,2])$dist
      wmat[i,1:ni] <- ifelse(dist<=ringdist,1,0)
    }
    wmat[upper.tri(wmat)] <- t(wmat)[upper.tri(wmat)]
  }

  if (kern=="rect")  { wgt <- function(psi) {.5 } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { .75*(1-psi^2) } }
  if (kern=="bisq")  { wgt <- function(psi) { (15/16)*((1-psi^2)^2) } }
  if (kern=="tcub")  { wgt <- function(psi) { (70/81)*((1 - abs(psi)^3)^3) } }
  if (kern=="trwt")  { wgt <- function(psi) { (35/32)*((1 - psi^2)^3) } }

  if (method=="kernel") {
    for (i in seq(1,n)) {
      dist <- geodistance(longvar=coormat[,1],latvar=coormat[,2],lotarget=coormat[i,1],latarget=coormat[i,2])$dist
      h <- quantile(dist,window)
      wmat[i,] <- ifelse(dist<=h,wgt(dist/h)/h,0)
    }
  }

  diag(wmat) = 0
  tmat <- rowSums(wmat)
  tmat <- ifelse(tmat==0,1,tmat)
  if (eigenvalues==TRUE) {
    dmat <- sqrt(1/tmat)
    dmat <- diag(dmat)%*%wmat%*%diag(dmat)
    eigvar <- eigen(dmat,only.values=TRUE,symmetric=TRUE)$values
  }
  wmat <- as.matrix(as.data.frame(wmat)/tmat)

  out <- list(wmat,eigvar)
  names(out) <- c("wmat","eigvar")
  return(out)
}


############################################################
# Functions processing radar images and their example. 
############################################################

######################################
# processing Radar images from V. Chen
######################################

## b727 <- matrix(scan("b727s.dat"),ncol=512,byrow=TRUE)
## b727 <- matrix(scan("b727r.dat"),ncol=256,byrow=TRUE)

## ncol <- dim(b727)[2]/2
## b727r<- b727[,-1+2*(1:ncol)]   
## b727i<- b727[,2*(1:ncol)]       

## B727 <- t(Conj(b727r + i*b727i))
## nrow <- dim(B727)[1]
## ncol <- dim(B727)[2]

## shift fft 1D data
## -----------------

fftshift <- function(x) {
  lng <- length(x)
  y <- numeric(lng)
  y[1:(lng/2)] <- x[(lng/2+1):(lng)]
  y[(lng/2+1):lng] <- x[1:(lng/2)]
  y
}

## The following tow S-routine concludes our Radar experiments
## -------------------------------------------------------------
## b727 <- matrix(scan("b727s.dat"),ncol=512,byrow=TRUE)
## b727 <- matrix(scan("b727r.dat"),ncol=256,byrow=TRUE)

## Generating the original image
## -----------------------------

## e.g. B727 <- showRadar(b727)

showRadar<- function(x, plot=TRUE) {
  ncol <- dim(x)[2]/2
  xr<- x[,-1+2*(1:ncol)]   
  xi<- x[,2*(1:ncol)]       
  y <- t(Conj(xr + 1i*xi))
  oy <- t(apply(apply(y,2,fftshift),1,fftshift))
  oy <- apply(oy,2,fft)
  oy <- t(apply(apply(oy,2,fftshift),1,fftshift))
  if(plot) image(Mod(t(oy)))
  oy
}

## Enhancing Radar Image by Gabor Transform (V.Chen)
## ------------------------------------------------

## gtime <- 128
## scale <- 50

## e.g. B727 <- cgtRadar(b727,gtime,scale)
## e.g. 
##      b727 <- matrix(scan("b727s.dat"),ncol=512,byrow=TRUE)
##      b727 <- matrix(scan("b727r.dat"),ncol=256,byrow=TRUE)
##      ncol <- dim(b727)[2]/2
##      b727r<- b727[,-1+2*(1:ncol)]   
##      b727i<- b727[,2*(1:ncol)]       

##      b727 <- t(Conj(complex(real=b727r, imag=b727i)))
##      B727 <- cgtRadar(b727,gtime,scale,flag=FALSE)
##      image(t(apply(Mod(B727$cgtout[,2,]),1,fftshift))) 

cgtRadar <- function(x,gtime,scale,flag=TRUE,plot=TRUE)
{
  y <- x
  if(flag) {
    ncol <- dim(x)[2]/2
    xr<- x[,-1+2*(1:ncol)]   
    xi<- x[,2*(1:ncol)]       
    y <- t(Conj(xr + 1i*xi))
  }
  nrow <- dim(y)[1]
  ncol <- dim(y)[2]
  cgty <- array(0+0i,c(ncol,nrow,gtime))
  cgtyy <- array(0+0i,c(ncol,nrow,gtime))
  for(k in 1:ncol) cgty[k,,] <-
    cgt(Re(y[,k]),gtime,2/gtime,scale,plot=FALSE) + 1i*cgt(Im(y[,k]),gtime,2/gtime,scale,plot=FALSE)
  for(k in 1:nrow)
    cgtyy[,k,] <- t(apply(cgty[,k,],1,fftshift))
  oy <- t(apply(apply(Mod(cgty),c(1,3),mean),1,fftshift))
  if(plot)
    image(oy)
  list(output=oy,cgtout=cgtyy)
}

## Generate X, Y, Z for the command : cloud

cloudXYZ <- function(dim1, dim2,dim3)
{
  II <- rep(1:dim1,dim2*dim3)
  dim(II) <- c(dim1,dim2,dim3)
  KK <- rep(1:dim3,dim2)
  dim(KK) <- c(dim3,dim2)
  KK <- t(KK)
  KK <- rep(KK,dim1)
  dim(KK) <- c(dim2 * dim3, dim1)
  KK <- t(KK)
  dim(KK) <- c(dim1, dim2, dim3)
  JJ <- rep(1:dim2,dim1)
  dim(JJ) <- c(dim2,dim1)
  JJ <- t(JJ)
  JJ <- rep(JJ,dim3)
  dim(JJ) <- c(dim1,dim2,dim3)
  list(II=II, JJ=JJ, KK=KK)
}

## Commands show points :

## cgtout <- B727$cgtout
## tt <-  cloudXYZ(dim(cgtout)[1], dim(cgtout)[2], dim(cgtout)[3])
## range(abs(cgtout))
## SS <- (abs(cgtout) > 4.0)
## XX <- tt$II[SS]
## YY <- tt$JJ[SS]
## ZZ <- tt$KK[SS]
## XX[length(XX)+1] <- 0
## XX[length(XX)+1] <- dim1
## YY[length(YY)+1] <- 0
## YY[length(YY)+1] <- dim2
## ZZ[length(ZZ)+1] <- 0
## ZZ[length(ZZ)+1] <- dim3
## cloud(ZZ~XX*YY)

## showpoints <- function(cgtout,low,high) {
## X <- 1:crossrange
## XX <- X %*% t(rep(1,crossrange))
## YY <- t(XX)
## freq <- dim(cgtout)[2]
## ZZ <- matrix(0,crossrange,crossrange)
## modmax <- max(Mod(cgtout))
## for(j in 1:crossrange)
##  for(k in 1:crossrange)
##    ZZ[j,k] <- maxpos(Mod(cgtout[j,,k]),low*modmax,high*modmax)
## cloud(ZZ~XX*YY)
##}

## maxpos <- function(x,low,high) {
##   pos <- 1       
##   lng <- length(x)
##   maxval <- max(x)
##   if((maxval >= low) & (maxval <= high)) { 
##    for(j in 1:lng) {
##      if(x[j] >= maxval) {
##        pos <- j
##      }
##    }
##   }
##  pos
## }

## showpoints <- function(cgtout,low,high) {
##  crossrange <- dim(cgtout)[1]
##  gabortime <- dim(cgtout)[3]
##  freq <- dim(cgtout)[2]
##  modmax <- max(Mod(cgtout))
##  tt <-persp(matrix(1,crossrange,gabortime),zlim=1:freq)
##  TT <- (Mod(cgtout) > low * modmax) * (Mod(cgtout) <= high * modmax)
##  for(k in 1:gabortime) {
##     for(j in 1:freq) {
##       for(i in 1:crossrange) { 
##         if(TT[i,j,k]==TRUE) 
##           points(perspp(i,k,j,tt))
##       }
##     }
##   }
##  tt
## }

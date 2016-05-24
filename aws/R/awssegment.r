#
#    R - function  aws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and  
#    Volatility models                                                         
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2006 Weierstrass-Institut fuer
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  Joerg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#
#     default parameters:  see function setawsdefaults
#       
aws.segment <- function(y, level,delta=0,hmax=NULL,hpre=NULL,
                varmodel="Constant",lkern="Triangle",scorr=0,ladjust=1,
                wghts=NULL,u=NULL,varprop=.1,ext=0,
                graph=FALSE,demo=FALSE,fov=NULL)
{
#
#    first check arguments and initialize
#
args <- match.call()
dy<-dim(y)
if(is.null(dy)) dy <- length(y)
mask <- rep(TRUE,length(y))
if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
if(!(varmodel %in% c("Constant","Linear","Quadratic"))) stop("Model for variance not implemented")
#
#   check for segmentation parameters
#
if(is.null(level)||is.null(delta)||level+delta>max(y)||level-delta<min(y)){
stop(paste("Inproper specifications for level ",level," or delta ",delta,
       "\n Values specified outside range of y"))
}
#
#   set appropriate defaults
#
if(is.null(wghts)) wghts <- c(1,1,1)
wghts <- switch(length(dy),c(0,0),c(wghts[1]/wghts[2],0),wghts[1]/wghts[2:3])
if(is.null(wghts)) wghts <- c(0,0)
cpar<-setawsdefaults(dy,mean(y),"Gaussian",lkern,"Uniform",TRUE,FALSE,ladjust,hmax,1,wghts)
lkern <- cpar$lkern
lambda <- 1.25*cpar$lambda # Gaussian case
maxvol <- cpar$maxvol
k <- cpar$k
kstar <- cpar$kstar
cpar$tau1 <- cpar$tau1*2 
cpar$tau2 <- cpar$tau2*2 
hmax <- cpar$hmax
shape <- cpar$shape
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- hmax*0.42445*4
    }
d <- cpar$d
n<-length(y)
#
#    set threshold
#
thresh <- setawsthresh(d,kstar,ladjust,ext)
beta <- switch(d,1.06,1.42,1.33)
#  optimised for sample sizes
#     d=1: n=c(2000,4000,8000)
#     d=2: n=c(256^2,512^2,1024^2)=c(65536,262144,1048576)
#     d=3: n=c(32^3,64^2*32,128^2*32)=c(32768,131072,524288)
# 
#   family dependent transformations 
#
zfamily <- awsgfamily(y,scorr,d)
sigma2 <- zfamily$sigma2
varest <- 1/sigma2
h0 <- zfamily$h0
rm(zfamily)
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- hmax*0.42445*4
    hinit <- 0.42445*4
    }
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
##  this is the version on a grid
n <- length(y)
n1 <- switch(d,n,dy[1],dy[1])
n2 <- switch(d,1,dy[2],dy[2])
n3 <- switch(d,1,1,dy[3])
if(is.null(fov)) fov <- n
interval <- if(delta>0) paste("Central interval: (",level-delta,",",level+delta,")") else paste("Level: ",level)
cat("Running segmentation algorithm with following parameters:\n",
interval,"\n",
"Critical value: ",thresh,"  Extension: ",ext,"  Field of view: ", fov,"\n")
#
#    Initialize  for the iteration
#
fix <- rep(FALSE,n)
zobj<-list(ai=y, bi0= rep(1,n), bi2=rep(1,n), bi=rep(1,n), theta=y/shape, fix=rep(FALSE,n)  )
segment <- array(0,c(n1,n2,n3))
vred<-rep(1,n)
mae<-NULL
lambda0<-1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
#
#   produce a presmoothed estimate to stabilize variance estimates
#
if(is.null(hpre)) hpre<-20^(1/d)
dlw<-(2*trunc(hpre/c(1,wghts))+1)[1:d]
hobj <- .Fortran("caws",as.double(y),
                       as.logical(fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.double(hpre),
                       as.double(rep(1,n)),
                       as.double(1e40),
                       as.double(zobj$theta),
                       bi=double(n),
                       double(n),
                       as.double(zobj$bi0),
                       ai=double(n),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(0.25),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","ai")]
hobj$theta <- hobj$ai/hobj$bi
dim(hobj$theta) <- dim(hobj$bi) <- dy
#
#   iteratate until maximal bandwidth is reached
#
total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
while (k<=kstar) {
      hakt0 <- gethani(1,10,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,10,lkern,1.25^k,wghts,1e-4)
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hakt <- hakt*0.42445*4
    }
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
if(scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
# Correction for spatial correlation depends on h^{(k)} 
hakt0<-hakt
# heteroskedastic Gaussian case
zobj <- .Fortran("segment",as.double(y),
                       fix=as.logical(fix),
                       as.double(level),
                       as.double(delta),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       bi=as.double(zobj$bi),
                       bi2=double(n),
                       bi0=as.double(zobj$bi0),
                       gi=double(n),
                       vred=double(n),
                       theta=as.double(zobj$theta),
                       as.integer(lkern),
                       as.double(0.25),
                       double(prod(dlw)),
                       as.double(wghts),
                       pvalues=double(n),# array for pvalues
                       as.integer(segment),# previous segmentation array 
                       segment=as.integer(segment),# new array for segment (-1,0,1)
                       as.double(beta),
                       as.double(thresh),
                       as.double(ext),
                       as.double(fov),
                       varest=as.double(varest),
                       PACKAGE="aws")[c("fix","bi","bi0","bi2","vred","pvalues",
                       "segment","theta","gi","hakt","varest")]
vred[!fix]<-zobj$vred[!fix]
if(hakt>n1/2) zobj$bi0 <- rep(max(zobj$bi),n)
pvalues <- zobj$pvalues
segment <- zobj$segment
varest <- zobj$varest
fix <- zobj$fix
dim(zobj$theta) <- dim(zobj$gi) <- dim(pvalues) <- dim(segment) <- dim(fix) <- dim(zobj$bi) <- dy
if(graph){
#
#     Display intermediate results if graph == TRUE
#
if(d==1){ 
oldpar<-par(mfrow=c(1,3),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y,ylim=range(y,zobj$theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(zobj$theta,lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(segment,type="l",main="Segmentation result",ylim=c(-1,1))
plot(zobj$bi,type="l",ylim=range(0,zobj$bi))
title("Sum of weights")
} 
if(d==2){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(zobj$theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(zobj$theta),3)," max=",signif(max(zobj$theta),3)))
image(zobj$gi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(zobj$gi),3)," mean=",signif(mean(zobj$gi),3)," max=",signif(max(zobj$gi),3)))
image(segment,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(-1,1))
title("Segmentation result")
}
if(d==3){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(zobj$theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(zobj$theta),3)," max=",signif(max(zobj$theta),3)))
image(zobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(zobj$bi),3)," mean=",signif(mean(zobj$bi),3)," max=",signif(max(zobj$bi),3)))
image(segment[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(-1,1))
title("Segmentation result")
} 
par(oldpar)
}
#
#    Calculate MAE and MSE if true parameters are given in u 
#    this is for demonstration and testing for propagation (parameter adjustments) 
#    only.
#
if(!is.null(u)) {
   cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    signif(mean((zobj$theta-u)^2),3),"   MAE: ",
		    signif(mean(abs(zobj$theta-u)),3)," mean(bi)=",
		    signif(mean(zobj$bi),3),"\n")
   mae<-c(mae,signif(mean(abs(zobj$theta-u)),3))
		    }
if(demo) readline("Press return")
#
#   Prepare for next iteration
#
#
#   Create new variance estimate
#
if(sum(zobj$fix)<prod(dy)/4){
# the estimates of sigma should be stable enough otherwise
# avoids estimating variances from small number of design points
vobj <- awsgsigma2(y,mask,hobj,zobj[c("theta","gi")],varmodel,varprop,h0)
sigma2 <- vobj$sigma2inv
coef <- vobj$coef
rm(vobj)
}
x<-1.25^(k-1)
scorrfactor<-x/(3^d*prod(scorr)*prod(h0)+x)
lambda0<-lambda*scorrfactor
if (max(total) >0) {
      cat("step:",k,"  hakt=",hakt,"  progress:",signif(total[k],2)*100,"% .  fixed=",sum(zobj$fix),"\n",sep="")
     }
k <- k+1
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
###   component var contains an estimate of Var(zobj$theta) 
###   
vartheta <- zobj$bi2/zobj$bi^2
vartheta<-vartheta/Spatialvar.gauss(hakt/0.42445/4,h0+1e-5,d)*Spatialvar.gauss(hakt/0.42445/4,1e-5,d)
awssegmentobj(y,zobj$theta,segment,vartheta,level,delta,hakt,1/sigma2,lkern,lambda,ladjust,TRUE,FALSE,
              args,homogen=FALSE,earlystop=FALSE,family="Gaussian",wghts=wghts,
              scorr=scorr,varmodel=varmodel,vcoef=coef,mae=mae)
}
setawsthresh <- function(d,kstar,ladjust,ext){
ladjust <- min(ladjust,switch(d,1.68,2.37,2.44))
#  for ladjust>=switch(d,1.68,2.37,2.44) use nonadaptive thresholds
switch(d,1.65  + 0.1952*log(kstar)-0.0473*ladjust-0.6771*ext,
         1.729 + 0.2831*log(kstar)-0.2715*ladjust-0.4576*ext,
         1.696 + 0.4010*log(kstar)-0.2975*ladjust-0.4273*ext)
}

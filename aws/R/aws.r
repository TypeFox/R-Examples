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
aws <- function(y,hmax=NULL,aws=TRUE,memory=FALSE,family="Gaussian",
                lkern="Triangle",homogen=TRUE,aggkern="Uniform",
                sigma2=NULL,shape=NULL,scorr=0,spmin=0.25,
                ladjust=1,wghts=NULL,u=NULL,graph=FALSE,demo=FALSE,
                testprop=FALSE,maxni=FALSE)
{
#
#   this version uses neighborhoods with an increase in potential 
#   variance reduction by a factor of 1.25 from one iteration step 
#   to the next
#
#    wghts is interpreted as voxel extensions ..., wghts for nonexisting dimensions are are set to INFTY
#
#    first check arguments and initialize
#
args <- match.call()
dy<-dim(y)
if(is.null(dy)) dy <- length(y)
if(length(dy)>3) stop("AWS for more than 3 dimensional grids is not implemented")
#
#   set appropriate defaults
#
if(is.null(wghts)) wghts <- c(1,1,1)
wghts <- switch(length(dy),c(0,0),c(wghts[1]/wghts[2],0),wghts[1]/wghts[2:3])
if(family=="NCchi"){
varstats <- sofmchi(shape/2) # precompute table of mean, sd and var for 
#
#   NCchi for noncentral chi with shape=degrees of freedom and theta =NCP
#
}
cpar<-setawsdefaults(dy,mean(y),family,lkern,aggkern,aws,memory,ladjust,hmax,shape,wghts)
lambda <- cpar$lambda
hmax <- cpar$hmax
shape <- cpar$shape
d <- cpar$d
n<-length(y)
# 
#   family dependent transformations that depend on the value of family
#
zfamily <- awsfamily(family,y,sigma2,shape,scorr,lambda,cpar)
cpar <- zfamily$cpar
lambda <- zfamily$lambda
sigma2 <- zfamily$sigma2
h0 <- zfamily$h0
y <- zfamily$y
lkern <- cpar$lkern
rm(zfamily)
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
##  this is the version on a grid
n1 <- switch(d,n,dy[1],dy[1])
n2 <- switch(d,1,dy[2],dy[2])
n3 <- switch(d,1,1,dy[3])
#
#    Initialize  for the iteration
#
maxvol <- cpar$maxvol
k <- cpar$k
kstar <- cpar$kstar
tobj<-list(bi= rep(1,n), bi2= rep(1,n), theta= y/shape, fix=rep(FALSE,n))
if(maxni) bi <- tobj$bi
zobj<-list(ai=y, bi0= rep(1,n))
hhom <- rep(1,n)
if(family=="Gaussian"&length(sigma2)==n) vred<-rep(1,n)
mae<-NULL
lambda0<-1e50 # that removes the stochstic term for the first step, Initialization by kernel estimates
if(testprop) {
#
#  prepare to  check for alpha in propagation condition (to adjust value of lambda using parameter ladjust)
#
       if(is.null(u)) u <- 0
       cpar <- c(cpar, list(n1=n1,n2=n2,n3=n3,n=n1*n2*n3,family=family,u=u))
       propagation <- NULL
    } 
#
#   iteratate until maximal bandwidth is reached
#
cat("Progress:")
total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
while (k<=kstar) {
      hakt0 <- gethani(1,1.25*hmax,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,1.25*hmax,lkern,1.25^k,wghts,1e-4)
      cat("step",k,"hakt",hakt,"\n")
if(lkern==5) {
#  assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hakt <- hakt*0.42445*4
    }
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
if(family=="Gaussian"&scorr[1]>=0.1) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,d)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,d)
# Correction for spatial correlation depends on h^{(k)} 
if(family=="Gaussian"&length(sigma2)==n){
# heteroskedastic Gaussian case
zobj <- .Fortran("chaws",as.double(y),
                       as.logical(tobj$fix),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi2=double(n),
                       bi0=double(n),
                       vred=double(n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","vred","ai","hakt")]
vred[!tobj$fix]<-zobj$vred[!tobj$fix]
} else {
# all other cases
if(cpar$mcode!=6){
   zobj <- .Fortran("caws",as.double(y),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       bi=as.double(tobj$bi),
                       bi2=double(n),
                       bi0=double(n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","ai","hakt","hhom")]
   } else {
   zobj <- .Fortran("caws6",as.double(y),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt=as.double(hakt),
                       hhom=as.double(hhom),
                       as.double(lambda0),
                       as.double(tobj$theta),
                       as.double(fncchiv(tobj$theta,varstats)/2),
                       bi=as.double(tobj$bi),
                       bi2=double(n),
                       bi0=double(n),
                       ai=as.double(zobj$ai),
                       as.integer(cpar$mcode),
                       as.integer(lkern),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       PACKAGE="aws",DUP=TRUE)[c("bi","bi0","bi2","ai","hakt","hhom")]
   }                     
}
if(family%in%c("Bernoulli","Poisson")) zobj<-regularize(zobj,family)
dim(zobj$ai)<-dy
tobj<-updtheta(zobj,tobj,cpar)
dim(tobj$theta)<-dy
if(maxni) bi <- tobj$bi <- pmax(bi,tobj$bi)
dim(tobj$bi)<-dy
dim(tobj$eta)<-dy
dim(tobj$fix)<-dy
if(homogen) hhom <- zobj$hhom
#
#  if testprop == TRUE
#  check alpha in propagation condition (to adjust value of lambda)
#  
if(testprop) propagation <- awstestprop(y,family,tobj,zobj,sigma2,hakt,cpar,u,propagation)
if(graph){
#
#     Display intermediate results if graph == TRUE
#
if(d==1){ 
oldpar<-par(mfrow=c(1,2),mar=c(3,3,3,.2),mgp=c(2,1,0))
plot(y,ylim=range(y,tobj$theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(tobj$theta,lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(tobj$bi,type="l",ylim=range(0,tobj$bi))
lines(tobj$eta*max(tobj$bi),col=2)
lines(hhom/max(hhom)*max(tobj$bi),col=3)
title("Sum of weights, eta and hhom")
} 
if(d==2){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(tobj$theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$fix,col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("Estimates fixed")
}
if(d==3){ 
oldpar<-par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
image(tobj$theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(tobj$theta),3)," max=",signif(max(tobj$theta),3)))
image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
image(tobj$fix[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n",zlim=c(0,1))
title("Estimates fixed")
} 
par(oldpar)
}
#
#    Calculate MAE and MSE if true parameters are given in u 
#    this is for demonstration and testing for propagation (parameter adjustments) 
#    only.
#
if(!is.null(u)) {
   cat("bandwidth: ",signif(hakt,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
                    signif(mean((tobj$theta-u)^2),3),"   MAE: ",
		    signif(mean(abs(tobj$theta-u)),3)," mean(bi)=",
		    signif(mean(tobj$bi),3),"mean hhom",signif(mean(hhom),3),"\n")
   mae<-c(mae,signif(mean(abs(tobj$theta-u)),3))
		    }
if(demo) readline("Press return")
#
#   Prepare for next iteration
#
x<-1.25^k
scorrfactor<-x/(3^d*prod(scorr)*prod(h0)+x)
lambda0<-lambda*scorrfactor
if (max(total) >0) {
      cat(signif(total[k],2)*100,"% . ",sep="")
     }
k <- k+1
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
###   component var contains an estimate of Var(tobj$theta) if aggkern="Uniform", or if !memory 
###   
if( family=="Gaussian"&length(sigma2)==n){
# heteroskedastic Gaussian case 
vartheta <- tobj$bi2/tobj$bi^2
#  pointwise variances are reflected in weights
} else {
vartheta <- switch(family,Gaussian=sigma2,
                          Bernoulli=tobj$theta*(1-tobj$theta),
                          Poisson=tobj$theta,
                          Exponential=tobj$theta^2,
                          Volatility=2*tobj$theta,
                          Variance=2*tobj$theta,0)*tobj$bi2/tobj$bi^2
vred<-tobj$bi2/tobj$bi^2
}
sigma2 <- switch(family,Gaussian=sigma2,
                        Bernoulli=tobj$theta*(1-tobj$theta),
                        Poisson=tobj$theta,
                        Exponential=tobj$theta^2,
                        Volatility=2*tobj$theta,
                        Variance=2*tobj$theta,0)
if( family=="Gaussian"){
vartheta<-vartheta/Spatialvar.gauss(hakt/0.42445/4,h0+1e-5,d)*Spatialvar.gauss(hakt/0.42445/4,1e-5,d)
}
awsobj(y,tobj$theta,vartheta,hakt,sigma2,lkern,lambda,ladjust,aws,memory,
   args,homogen,earlystop=FALSE,family=family,wghts=wghts,mae=mae,ni=tobj$bi)
}
#######################################################################################
#
#        Auxilary functions
#
#######################################################################################
#
#        Set default values
#
# default values for lambda and tau are chosen by propagation condition
# (strong version) with alpha=0.05 (Gaussian) and alpha=0.01 (other family models)
# see script aws_propagation.r
#
#######################################################################################
setawsdefaults <- function(dy,meany,family,lkern,aggkern,aws,memory,ladjust,hmax,shape,wghts){
lkern<-switch(lkern,Triangle=2,Quadratic=3,Cubic=4,Plateau=1,Gaussian=5,2)
#
#   univariate case
#
if(is.null(dy)){
      d<-1
} else {
      d <- length(dy)
}
if(is.null(hmax)) hmax <- switch(d,250,12,5)
if(aws) lambda <- ladjust*switch(family,
# old values                  Gaussian=switch(d,11.3,6.1,6.2),# see inst/scripts/adjust.r for alpha values
Gaussian=switch(d,8.1,5.4,4.9),# see ladjgaussx.r in R/aws/ladj
# old values          Bernoulli=switch(d,9.6,7.6,6.9),
Bernoulli=switch(d,5.3,5,4.8),
                  Exponential=switch(d,7.1,6.1,5.5),
                  Poisson=switch(d,7.7,6.1,5.9),
                  Volatility=switch(d,5,4,3.7),
                  Variance=switch(d,12.8,6.1,6.1),
                  NCchi=switch(d,16.8,16.8,16.8),# worst case theta=3
                  switch(d,11.3,6.1,.96)) else lambda <- 1e50
#
#    determine heta for memory step
#
heta <- switch(family,   Gaussian=1.25,
               Bernoulli=5/meany/(1-meany),
               Exponential=20,
               Poisson=max(1.25,20/meany),
               Volatility=20,
               Variance=20,
               NCchi=1.25,1.25)^(1/d)
ktau <- log(switch(d,100,15,5))
#
# stagewise aggregation 
#
if(!aws){
#
# force aggkern = "Triangle"  here
#
#  we need a tau for stagewise aggregation that fulfils the propagation condition
#
    aggkern <- "Triangle"
    cat("Stagewise aggregation: Triangular aggregation kernel is used\n")
#
#   this is the first bandwidth to consider for stagewise aggregation
#
    if(!memory){
        heta <- hmax
        cat("Neither PS nor Stagewise Aggregation is specified, compute kernel estimate with bandwidth hmax\n") 
    } else {
    tau <- switch(family,
                  Gaussian=switch(d,.28,.21,.21),
                  Bernoulli=switch(d,.36,.36,.36),
                  Exponential=switch(d,1.3,1.3,1.3),
                  Poisson=switch(d,.36,.21,.21),
                  Volatility=switch(d,1.3,1.3,1.3),
                  Variance=switch(d,1.3,1.3,1.3),
                  NCchi=switch(d,.28,.21,.21),
                  switch(d,.28,.21,.21))
    }
} else {
#
#   set appropriate value for tau (works for all families)
#
if(memory) {
tau <- 8 
} else {
tau <- 1e10
heta <- 1e10
}
}
if(memory) tau1<-ladjust*tau else tau1 <- 1e10
#
#    adjust for different aggregation kernels
#
if(aggkern=="Triangle") tau1<-2.5*tau1
tau2<-tau1/2
#
#   set maximal bandwidth
#
# uses a maximum of about 500, 450 and 520  points, respectively.
mcode <- switch(family,
                Gaussian=1,
                Bernoulli=2,
                Poisson=3,
                Exponential=4,
                Volatility=4,
                Variance=5,
                NCchi=6,-1)
if(mcode < 0) stop(paste("specified family ",family," not yet implemented"))
lambda<-lambda*1.8
if(is.null(shape)) shape<-1
maxvol <- getvofh(hmax,lkern,wghts)
kstar <- as.integer(log(maxvol)/log(1.25))
if(aws||memory) k <- switch(d,1,3,6) else k <- kstar
if(aws) cat("Running PS with lambda=",signif(lambda,3)," hmax=",hmax,"number of iterations:",kstar," memory step",if(memory) "ON" else "OFF","\n")
else cat("Stagewise aggregation \n")
list(heta=heta,tau1=tau1,tau2=tau2,lambda=lambda,hmax=hmax,d=d,mcode=mcode,shape=shape,
     aggkern=aggkern,ktau=ktau,kstar=kstar,maxvol=maxvol,k=k,lkern=lkern,wghts=wghts)
}
#######################################################################################
#
#    IQRdiff (for robust variance estimates
#
#######################################################################################
IQRdiff <- function(y) IQR(diff(y))/1.908
#######################################################################################
#
#    Kullback-Leibler distances
#
#######################################################################################
KLdist <- function(mcode,th1,th2,bi0,shape){
   th12<-(1-0.5/bi0)*th2+0.5/bi0*th1
   z<-switch(mcode,(th1-th2)^2,
                th1*log(th1/th12)+(1.-th1)*log((1.-th1)/(1.-th12)),
                th1*log(th1/th12)-th1+th12,
                th1/th2-1.-log(th1/th2),
                shape/2*(th2/th1-1.)+(shape/2-1)*log(th1/th2))
   z[is.na(z)]<-0
   z
}
KLdist6 <- function(th1,th2,bi0,shape1,shape2){
# case mcode=6 Symmetrized Gauss approximation
(th1-th2)^2/(shape1+shape2)
   z[is.na(z)]<-0
   z
}
####################################################################################
#
#    Memory step for local constant aws
#
####################################################################################
updtheta<-function(zobj,tobj,cpar){
heta<-cpar$heta
hakt<-zobj$hakt
bi<-zobj$bi
bi2<-zobj$bi2
thetanew<-zobj$ai/bi
if(hakt>heta) {
#
#   memory step
#
mcode<-cpar$mcode
shape<- cpar$shape
aggkern <- cpar$aggkern
tau1<-cpar$tau1
tau2<-cpar$tau2
kstar<-cpar$ktau
tau<-2*(tau1+tau2*max(kstar-log(hakt),0))
theta<-tobj$theta
thetanew[tobj$fix]<-theta[tobj$fix]
if(mcode<6){
eta<-switch(aggkern,
            "Uniform"=as.numeric(zobj$bi0/tau*
              KLdist(mcode,thetanew,theta,max(zobj$bi0),shape)>1),
            "Triangle"=pmin(1,zobj$bi0/tau*
              KLdist(mcode,thetanew,theta,max(zobj$bi0),shape)),
             as.numeric(zobj$bi0/tau*
              KLdist(mcode,thetanew,theta,max(zobj$bi0),shape)>1))
} else {
  eta <- 1
#
#  no memory step implemented in this case
#
}
eta[tobj$fix]<-1
bi <- (1-eta)*bi + eta * tobj$bi
bi2 <- (1-eta)*bi2 + eta * tobj$bi2
thetanew <- (1-eta)*thetanew + eta * theta
} else {
#
#  no memory step
#
eta <- rep(0,length(thetanew))
}
list(theta=thetanew,bi=bi,bi2=bi2,eta=eta,fix=(tobj$fix|eta==1))
}
####################################################################################
#
#    Regularize for Bernoulli and Poisson models if parameter estimates are
#    at the border of their domain
#
####################################################################################
regularize <- function(zobj,family){
if(family%in%c("Bernoulli")){
   zobj$ai <- .1/zobj$bi+zobj$ai
   zobj$bi <- .2/zobj$bi+zobj$bi
}
if(family%in%c("Poisson")) zobj$ai <- 0.1/zobj$bi+zobj$ai
zobj
} 
############################################################################
#
#   transformations that depend on the specified family
#
############################################################################
awsfamily <- function(family,y,sigma2,shape,scorr,lambda,cpar){
h0 <- 0
if(family=="Gaussian") {
  d <- cpar$d
  if(any(scorr>0)) {
         h0<-numeric(length(scorr))
         for(i in 1:length(h0))
         h0[i]<-geth.gauss(scorr[i])
         if(length(h0)<d) h0<-rep(h0[1],d)
         cat("Corresponding bandwiths for specified correlation:",h0,"\n")
}
    if(is.null(sigma2)) {
        sigma2 <- IQRdiff(as.vector(y))^2
        if(scorr[1]>0) sigma2<-sigma2*Varcor.gauss(h0)
        cat("Estimated variance: ", signif(sigma2,4),"\n")
    }
    if(length(sigma2)==1){
#   homoskedastic Gaussian case
    lambda <- lambda*sigma2*2 
    cpar$tau1 <- cpar$tau1*sigma2*2 
    cpar$tau2 <- cpar$tau2*sigma2*2 
    } else {
#   heteroskedastic Gaussian case
    if(length(sigma2)!=length(y)) 
       stop("sigma2 does not have length 1 or same length as y")
    lambda <- lambda*2 
    cpar$tau1 <- cpar$tau1*2 
    cpar$tau2 <- cpar$tau2*2 
    sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
    }
}
#
#   specify which statistics are needed and transform data if necessary
#
if(family=="Volatility"){
   family <- "Exponential"
   y <- y^2
   lambda <- 2*lambda 
# this accounts for the additional 1/2 in Q(\hat{theta},theta)
}
if(family=="Variance"){
   lambda <- 2*lambda/cpar$shape 
   cpar$tau1 <- cpar$tau1*2/cpar$shape 
   cpar$tau2 <- cpar$tau2*2/cpar$shape 
}
list(cpar=cpar,lambda=lambda,y=y,sigma2=sigma2,h0=h0)
}

getvofh <- function(bw,lkern,wght){
.Fortran("getvofh",
         as.double(bw),
         as.integer(lkern),
         as.double(wght),
         vol=double(1),
         PACKAGE="aws")$vol
}
gethani <- function(x,y,lkern,value,wght,eps=1e-2){
.Fortran("gethani",
         as.double(x),
         as.double(y),
         as.integer(lkern),
         as.double(value),
         as.double(wght),
         as.double(eps),
         bw=double(1),
         PACKAGE="aws")$bw
}

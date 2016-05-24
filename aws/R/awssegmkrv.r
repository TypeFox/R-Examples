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
aws.segment.krv <- function(y,hmax=NULL,sigma2=NULL,ladjust=1,wghts=NULL){
#
#    first check arguments and initialize
#
args <- match.call()
dy<-dim(y)
if(is.null(dy)){
      d<-1
} else {
      d <- length(dy)
}
if(d>3) stop("AWS for more than 3 dimensional grids is not implemented")
#
#   set appropriate defaults
#
if(is.null(wghts)) wghts <- c(1,1,1)
wghts <- switch(d,c(0,0),c(wghts[1]/wghts[2],0),wghts[1]/wghts[2:3])
if(is.null(hmax)) hmax <- switch(d,250,12,5)
maxvol <- getvofh(hmax,1,wghts)
lambda <- ladjust*switch(d,19,19,21)
kstar <- as.integer(log(maxvol)/log(1.25))
fov <- length(y)
#
#    set threshold
#
n <- length(y)
n1 <- switch(d,n,dy[1],dy[1])
n2 <- switch(d,1,dy[2],dy[2])
n3 <- switch(d,1,1,dy[3])
if(is.null(sigma2)) {
        sigma2 <- IQRdiff(as.vector(y))^2
	cat("Estimated variance: ", signif(sigma2,4),"\n")
	}
if(length(sigma2)!=n) sigma2 <- rep(sigma2[1],n)
varest <- 1/sigma2
#
#    Initialize  for the iteration
#
zobj<-list(bi=varest, th=y)
lambda0<-1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
#
#   produce a presmoothed estimate to stabilize variance estimates
#
#
#   iteratate until maximal bandwidth is reached
#
total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
k <- 1
krval <- matrix(0,2,kstar)
mae <- numeric(kstar)
while (k<=kstar) {
      hakt0 <- gethani(1,10,1,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,10,1,1.25^k,wghts,1e-4)
dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:d]
hakt0<-hakt
# heteroskedastic Gaussian case
zobj <- .Fortran("segmgkrv",as.double(y),
                       as.double(varest),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$th),
                       th=double(n),
                       bi=as.double(zobj$bi),
		       double(prod(dlw)),
		       as.double(wghts),
                       as.double(fov),
                       maxvalue=double(1),
                       minvalue=double(1),
		       PACKAGE="aws",DUP=FALSE)[c("th","bi","minvalue","maxvalue")]
dim(zobj$th) <- dim(zobj$bi) <- dy
krval[1,k] <- -zobj$minvalue
krval[2,k] <- zobj$maxvalue
#
#    Calculate MAE and MSE if true parameters are given in u 
#    this is for demonstration and testing for propagation (parameter adjustments) 
#    only.
#
mae[k] <- mean(abs(zobj$th))
#
#   Prepare for next iteration
#
#
#   Create new variance estimate
#
lambda0<-lambda
k <- k+1
gc()
}
cat("\n")
###                                                                       
###            end iterations now prepare results                                                  
###                                 
###   component var contains an estimate of Var(zobj$theta) 
###   
list(th=zobj$th,bi=zobj$bi,krval=krval,mae=mae)
}

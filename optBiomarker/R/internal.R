##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License 
## as published by the Free Software Foundation; either version 2 
## of the License, or (at your option) any later version.
##  
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
## Public License for more details.
##   
## You may also obtain a copy of the GNU General Public License from 
## the Free Software Foundation by visiting their Web site or by writing to
##   
##
##  Free Software Foundation, Inc.
##  59 Temple Place - Suite 330
##  Boston, MA 02111-1307
##  USA
##
################################################################################


################################################################################
##                          *** tnormAvg ***                        
################################################################################

## Function for calculating the expected value of  a truncated normal distribution
## -------------------------------------------------------------------------------

## Details:

##The truncated normal distribution has density:

## Suppose that Y~N(mu, sigma^2), and A=[a1,a2] is a subset of -Inf<y<Inf. 
## Now the conditional distribution of Y given A is called truncated normal
## distribution.

##    f(y|A) = (1/sigma) phi((y-mu)/sigma) / (Phi((a2-mu)/sigma) - Phi((a1-mu)/sigma))
##      for a1 <= y <= a2, and 0 otherwise.

##  mu is the mean of the original Normal distribution before truncation,
##  sigma is the corresponding standard deviation,
## a2 is the upper truncation point,
## a1 is the lower truncation point,
## phi(x) is the density of the standard normal distribution, and
## Phi(x) is the distribution function of the standard normal distribution. 

## The expected value of this truncated distribution is give by:

##  E[Y|A]=mu-sigma\fraction{phi((a2-mu)/sigma)-phi((a1-mu)/sigma))}{(Phi((a2-mu)/sigma) - Phi((a1-mu)/sigma))} 
##  For given value of E[Y|A], assuming that mu and sigma are known, the above is an equation with two unknowns
##  a1 and a2, and therefore can not be solved uniquely.

## This function computes this expectation as  function of a1, a2, mu and sigma


## Function for computing E[Y|A]
## ---------------------------------


tnormAvg<-function(a1=-Inf, a2=Inf,mu=0,sigma=1)
{    
    E<-mu-sigma*((dnorm((a2-mu)/sigma) - dnorm((a1-mu)/sigma))/(pnorm((a2-mu)/sigma)-pnorm((a1-mu)/sigma)))
    return(E)    
}



################################################################################
##                          *** objFun ***                        
################################################################################

## This function defines the objective function for obtaining the lower truncation
## point of a truncated normal distribution, given the average of the truncated
## distribution
## -------------------------------------------------------------------------------

## Details:

##The truncated normal distribution has density:

## Suppose that Y~N(mu, sigma^2), and A=[a1,a2] is a subset of -Inf<y<Inf. 
## Now the conditional distribution of Y given A is called truncated normal
## distribution.

##    f(y|A) = (1/sigma) phi((y-mu)/sigma) / (Phi((a2-mu)/sigma) - Phi((a1-mu)/sigma))
##      for a1 <= y <= a2, and 0 otherwise.

##  mu is the mean of the original Normal distribution before truncation,
##  sigma is the corresponding standard deviation,
## a2 is the upper truncation point,
## a1 is the lower truncation point,
## phi(x) is the density of the standard normal distribution, and
## Phi(x) is the distribution function of the standard normal distribution. 

## The expected value of this truncated distribution is give by:

##  E[Y|A]=mu-sigma\fraction{phi((a2-mu)/sigma)-phi((a1-mu)/sigma))}{(Phi((a2-mu)/sigma) - Phi((a1-mu)/sigma))} 

## The objective function is defined as (tnormAvg(para)-truncAvg)^2 to be minimized over para

objFun<-function(para, mu=0,sigma=1,truncAvg=log2(2))
{    
    (tnormAvg(a1=para,a2=Inf, mu=mu, sigma=sigma)-truncAvg)^2
}


################################################################################
##                          *** avg2Lower ***                        
################################################################################

## Function for calculating the lower truncation point of a normal distribution 
## for a given average of the truncated distribution with known values of mu
## and sigma
## ----------------------------------------------------------------------------

## Details:

##The truncated normal distribution has density:

## Suppose that Y~N(mu, sigma^2), and A=[a1,a2] is a subset of -Inf<y<Inf. 
## Now the conditional distribution of Y given A is called truncated normal
## distribution.

##    f(y|A) = (1/sigma) phi((y-mu)/sigma) / (Phi((a2-mu)/sigma) - Phi((a1-mu)/sigma))
##      for a1 <= y <= a2, and 0 otherwise.

##  mu is the mean of the original Normal distribution before truncation,
##  sigma is the corresponding standard deviation,
## a2 is the upper truncation point,
## a1 is the lower truncation point,
## phi(x) is the density of the standard normal distribution, and
## Phi(x) is the distribution function of the standard normal distribution. 

## The expected value of this truncated distribution is give by:

##  E[Y|A]=mu-sigma\fraction{phi((a2-mu)/sigma)-phi((a1-mu)/sigma))}{(Phi((a2-mu)/sigma) - Phi((a1-mu)/sigma))} 
##  For given value of E[Y|A], assuming that mu and sigma are known, the above is an equation with two unknowns
##  a1 and a2, and therefore can not be solved uniquely.

## This function solves the equation numerically for a1 assuming a2->Inf.

avg2Lower<-function(truncAvg=log2(2), mu=0,sigma=1)
{
    if(truncAvg <0 || truncAvg > 8)
        stop("truncAvg outside [0,8] not allowed")
    if(mu<0)stop("negative mu is not allowed")
    result<-optimize(f=objFun,interval=c(-10,10), mu=mu, sigma=sigma,truncAvg=truncAvg)
    lower<-result$minimum
    fmin<-result$objective
    if(is.finite(fmin)&& abs(fmin)>0.001)
        warning("value of objective is not near zero at solution")
    return(lower)
}




################################################################################
##               *** classpredict.lda, classpredict.knn ***                        
################################################################################

## Force predict to return class labels only (LDA and KNN)
## RF and SVM do it by default

classpredict.lda <- function(object, newdata)
predict(object, newdata = newdata)$class

classpredict.knn <- function(object, newdata) 
                   predict.ipredknn(object, newdata, type="class")


################################################################################
##                            *** dimSelect ***                        
################################################################################

## Function for selections of dimensions for interpolation

dimSelect<-function(x,val)
{

extr<-range(x)

if( val <extr[1] || val>extr[2]) stop("val must be within the range of x")

y<-abs(x-val)
id1<-order(y)[1]
if(x[id1]!=val){
if (x[id1]<val) idx<-c(id1,id1+1) else idx<-c(id1-1,id1)
} else {
idx<-c(id1,id1)
}
idx
}


################################################################################
##                            *** plot3dFun ***                        
################################################################################

## Function for ploting 3D-graph

plot3dFun<-function(panel,...)
{ ## begin function
with(panel,
{ ##with block

## Indices of database

## Produce the message that LDA is not yet implemented
if(method=="LDA"){ ## Message for LDA
rp.messagebox("LDA will be implemented in the next release of the package", title = "Message")
return()
} ## end ## Message for LDA



method.val<-c("RF", "SVM", "KNN")
idx.method<-(1:length(method.val))[method.val %in% method]

nbiom.val<-c(1:5,7,9,11,15,20,30,40,50,100)
nbiom.ex<-c(1:100)

nTrain.val<-c(10,20,50,100,250)
idx.nTrain<-dimSelect(nTrain.val, nTrain)

##idx.nTrain<-order(abs(nTrain.val-nTrain))[1]
nTrain.ex<-c(10,20,30,40,50,60,70,80,90,100,115, 130,145,160,175,190, 205,220,235,250)



##idx.nTrain<-order(abs(nTrain.ex-nTrain))[1]
i.nTrain<-order(abs(nTrain.ex-nTrain))[1]


sdB.val<-c(0.5,1.0,1.5,2.5)
idx.sdB<-dimSelect(sdB.val, sdB)
#idx.sdB<-order(abs(sdB.val-sdB))[1]
sdB.ex<-signif(seq(0.5,2.5,length=20),digits=2)
i.sdB<-order(abs(sdB.ex-sdB))[1]


sdW.val<-c(0.1,0.5,1.0,1.5)
idx.sdW<-dimSelect(sdW.val, sdW)
##idx.sdW<-order(abs(sdW.val-sdW))[1]
sdW.ex<-signif(seq(0.1,1.5,length=20),digits=2)
##idx.sdW<-order(abs(sdW.ex-sdW))[1]


foldAvg.val<-c(1.74,2.88,4.03,6.33)
idx.foldAvg<-dimSelect(foldAvg.val, foldAvg)
#idx.foldAvg<-order(abs(foldAvg.val-foldAvg))[1]
foldAvg.ex<-signif(seq(1.74,6.33,length=20),digits=3)
##idx.foldAvg<-order(abs(foldAvg.ex-foldAvg))[1]

nRep.val<-c(1,3,5,7,10)
idx.nRep<-dimSelect(nRep.val,nRep)
#idx.nRep<-order(abs(nRep.val-nRep))[1]
nRep.ex<-(1:10)
##idx.nRep<-order(abs(nRep.ex-nRep))[1]


##rho.val<-signif(seq(0, 0.95, length=20), digits=3)
##idx.rho<-order(abs(rho.val-rho))[1]


x<-seq(0,1,length=100)
y<-seq(0,0.5,length=20)
xat<-x[c(10,seq(20,100,by=20))]
xlab<-as.character(c(10,seq(20,100,by=20)))
yat<-y[c(1,seq(5,20,by=5))]
ylab.sdB<-as.character(sdB.ex[c(1,seq(5,20,by=5))])


##if (yaxis=="Biological variation"){ ## if yaxis=="Biological variation"


ydata<-error[,idx.nRep,,idx.sdW,idx.foldAvg,idx.nTrain,idx.method]

## Interpolate y for nTrain

if(diff(nTrain.val[idx.nTrain])!=0){
fc.nTrain<-(nTrain-nTrain.val[idx.nTrain[1]])/diff(nTrain.val[idx.nTrain])
ydata1<-ydata[,,,,,1]
ydata2<-ydata[,,,,,2]
ydata<-ydata1+fc.nTrain*(ydata2-ydata1)
} else {
ydata<-ydata[,,,,,1]
}

## Interpolate new y for foldAvg

if(diff(foldAvg.val[idx.foldAvg])!=0){
fc.foldAvg<-(foldAvg-foldAvg.val[idx.foldAvg[1]])/diff(foldAvg.val[idx.foldAvg])
ydata1<-ydata[,,,,1]
ydata2<-ydata[,,,,2]
ydata<-ydata1+fc.foldAvg*(ydata2-ydata1)
} else{
ydata<-ydata[,,,,1]
}

## Interpolate new y for sdW

if(diff(sdW.val[idx.sdW])!=0) {
fc.sdW<-(sdW-sdW.val[idx.sdW[1]])/diff(sdW.val[idx.sdW])
ydata1<-ydata[,,,1]
ydata2<-ydata[,,,2]
ydata<-ydata1+fc.sdW*(ydata2-ydata1)
} else {
ydata<-ydata[,,,1]
}

## Interpolate new y for repl
if(diff(nRep.val[idx.nRep])!=0){
fc.nRep<-(nRep-nRep.val[idx.nRep[1]])/diff(nRep.val[idx.nRep])
ydata1<-ydata[,1,]
ydata2<-ydata[,2,]
ydata<-ydata1+fc.nRep*(ydata2-ydata1)
} else {
ydata<-ydata[,1,]
}

ydata<-matapprox(ydata, rx=nbiom.val, cx=sdB.val, rout=nbiom.ex, cout=sdB.ex)

rglids<-rgl.ids()[,1]
n1ids<-length(rglids)
n.keep<-6
n.rmv<-12

if (n1ids>n.keep) rgl.pop(id=rglids[(n.keep+1):n1ids])
if(nids<-length(rgl.ids()[,1])==0){ ## start if (rgl.ids()[,1])==0)
p1<-surface3d(x,y,ydata, front="lines", col="green", size=1, alpha=0)
bbox3d(xat=xat,xlab=xlab,
         yat=yat,ylab=ylab.sdB,
         color=c("#333377","black"), emission="#333377", 
          specular="#3333FF", shininess=128, alpha=0.8)
title3d("","","No. of Biomarker","Biological variation","Error rate", cex=1.25)
}## end if (rgl.ids()[,1])==0)
 
p2<-surface3d(x,y,matrix(errorTol,ncol=20,nrow=100), col="red", alpha=0.2)
p3<-surface3d(x,y,ydata, front="lines", col="green", size=1, alpha=0.8)

zi<-ydata[,i.sdB]
test<-zi<=errorTol 
if (any(test)){
zopt<-(zi[test])[1]
p.opt<-(1:length(zi))[zi==zopt][1]
opt.num<-p.opt
} else {
p.opt<-NULL
opt.num<-"Not found"
}
if (!is.null(p.opt)){ ## start if (!is.null(p.opt))
lines3d(x[p.opt],y,ydata[p.opt,],color="red",front="lines",size=5,alpha=0.4)
lines3d(x,y[i.sdB],ydata[,i.sdB],color="red",front="lines",size=5,alpha=0.4)
p4<-spheres3d(x[p.opt],y[i.sdB],ydata[p.opt,i.sdB], color="red",radius=1/30,alpha=0.5)
p5<-text3d(x[p.opt],y[i.sdB],ydata[p.opt,i.sdB], paste(p.opt), cex=2,color="blue", alpha=0.8)
} ## end if (!is.null(p.opt))
title3d(paste("Optimal number:",opt.num),"","","","", cex=1.25)

##} ## end if yaxis=="Biological variation"

} #end of with block
) #end of with
panel
} ## end function


################################################################################
##                            *** yapprox ***                        
################################################################################

## Modify the approx function--

## This is a modified version of the function 'approx' (stats) 
## 1. The positions of x and y in the argument list has been
##    interchanged (to be able to use in apply).
## 2. The output line is modified to give 'y' values only
## 

yapprox<-function(y, x, xout, method="linear", n=50,
          yleft, yright, rule = 1, f = 0, ties = mean)
{
	result<-approx(x=x, y = y, xout=xout, method=method, n=n,
          yleft=yleft, yright=yright, rule = rule, f = f, ties = ties)$y
return(result)
}

################################################################################
##                            *** matapprox ***
################################################################################

## Define matrix version of the 'approx' function
## Function for interpolating data matrix

matapprox<-function(mat, rx, cx, nr=100,nc=20,rout, cout)
{

if (missing(rx))rx<-1:NROW(mat)
if (missing(cx))cx<-1:NCOL(mat)
nrx<-length(rx)
ncx<-length(cx)

if (missing(rout)) {
        if (nr<= 0)             stop("'matapprox' requires nr >= 1")
        rout <- seq(rx[1], rx[nrx], length.out= nr)
    }

if (missing(cout)) {
        if (nc<= 0)
            stop("'matapprox' requires nc >= 1")
        cout <- seq(cx[1], cx[ncx], length.out= nc)
    }

rmat<-apply(mat,2,yapprox,x=rx,xout=rout)
crmat<-apply(rmat,1,yapprox,x=cx,xout=cout)
return(t(crmat))
}

################################################################################
## End of:                         internal.R
################################################################################



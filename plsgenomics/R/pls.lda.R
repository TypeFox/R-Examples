### pls.lda.R  (2005-04-06)
###
###     Classification with PLS Dimension Reduction and Linear Discriminan Analysis
###
### Copyright 2004-04 Anne-Laure Boulesteix and Korbinian Strimmer
###
### 
###
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

pls.lda<-function(Xtrain, Ytrain, Xtest=NULL, ncomp, nruncv=0, alpha=2/3, priors=NULL)
{
ntrain<-nrow(Xtrain)
Ytrain<-as.factor(Ytrain)

if (is.vector(Xtest))
 {
 Xtest<-matrix(Xtest,1,length(Xtest))
 }
if (is.null(Xtest))
 {
 Xtest<-Xtrain
 }
if (nruncv==0&length(ncomp)>1) 
 stop("Since length(ncomp)>1, nruncv must be >0")
 
if (nruncv>0)
 {
 ncomp<-pls.lda.cv(Xtrain,Ytrain,ncomp=ncomp,nruncv=nruncv,alpha=alpha,priors=priors)
 }

pls.out<-pls.regression(Xtrain=Xtrain,Ytrain=transformy(Ytrain),Xtest=NULL,ncomp=ncomp)

Ztrain<-as.data.frame(matrix(pls.out$T,ntrain,ncomp))
Ztrain$y<-Ytrain
Ztest<-as.data.frame(scale(Xtest,center=pls.out$meanX,scale=FALSE)%*%pls.out$R)
if (is.null(priors))
     {
     lda.out <- lda(formula = y ~ ., data = Ztrain)
     }
    else
     {
     lda.out <- lda(formula = y ~ ., data = Ztrain, prior = priors)
     }


predclass<-predict(object=lda.out,newdata=Ztest)$class

return(list(predclass=predclass,ncomp=ncomp))
}

############################

transformy<-function(y)
{
y<-as.numeric(y)
K<-max(y)
if (K>2)
 {
 Y<-matrix(0,length(y),K)
 for (k in 1:K)
  {
  Y[,k]<-as.numeric(y==k)
  Y[,k]<-Y[,k]-mean(Y[,k])
  }
 }
else
 {
 Y<-matrix(y-mean(y),length(y),1)
 }
 
Y
}




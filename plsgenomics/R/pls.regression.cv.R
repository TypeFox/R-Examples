### pls.regression.cv.R  (2005-04-06)
###
###    Determination of the number of latent components to be used in PLS regression
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



pls.regression.cv<-function(Xtrain,Ytrain,ncomp,nruncv=20,alpha=2/3)
{
n<-nrow(Xtrain)
ntrain<-floor(n*alpha)
if (is.vector(Ytrain))
 {
 Ytrain<-matrix(Ytrain,n,1)
 }
samp<-matrix(0,ntrain,nruncv)


if (length(ncomp)==1)
 {
 if (ncomp==1)
  return(ncomp)
 else 
 ncomp<-1:ncomp
 }
  
for (i in 1:nruncv)
 {
 samp[,i]<-sample(n,ntrain)
 }
samp<-as.data.frame(samp)

errorcv<-sapply(samp,FUN=pls.regression.sample,Xtrain,Ytrain,ncomp)
meanerror<-apply(errorcv,MARGIN=1,FUN=mean)

ncomp<-ncomp[which.min(meanerror)]

return(ncomp)
}


#####################################

pls.regression.sample<-function(samp,X,Y,ncomp)
{
lncomp<-length(ncomp)
errorcv<-numeric(lncomp)
pls.out<-pls.regression(X[samp,],Y[samp,],ncomp=ncomp,Xtest=X[-samp,])

for (j in 1:lncomp)
  {
  errorcv[j]<-sum((pls.out$Ypred[,,j]-Y[-samp,])^2)
  }

errorcv
}

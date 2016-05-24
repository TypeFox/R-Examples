### pls.lda.cv.R  (2005-04-06)
###
###     Determination of the number of latent components to be used for Classification with PLS Dimension Reduction and Linear Discriminant Analysis
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



pls.lda.cv<-function(Xtrain,Ytrain,ncomp,nruncv=20,alpha=2/3,priors=NULL)
{

if (length(ncomp)==1)
 {
 if (ncomp==1)
  return(ncomp)
 else
  ncomp<-1:ncomp
 }
n<-nrow(Xtrain)
ntrain<-floor(n*alpha)
samp<-matrix(0,ntrain,nruncv)
for (i in 1:nruncv)
 {
 samp[,i]<-sample(n,ntrain)
 }
samp<-as.data.frame(samp)
errorcv<-sapply(samp,FUN=pls.lda.sample,Xtrain,Ytrain,ncomp=ncomp,priors=priors)
meanerror<-apply(errorcv,MARGIN=1,FUN=mean)
ncomp<-ncomp[which.min(meanerror)]
return(ncomp)
}


#####################################

pls.lda.sample<-function(samp,X,Y,ncomp,priors=NULL)
{
errorcv<-numeric(length(ncomp))

for (j in ncomp)
  {
  pls.lda.out<-pls.lda(Xtrain=X[samp,],Ytrain=Y[samp],Xtest=X[-samp,],ncomp=j,nruncv=0,priors=priors)
  errorcv[j]<-sum(pls.lda.out$predclass!=Y[-samp])
  }
errorcv
}

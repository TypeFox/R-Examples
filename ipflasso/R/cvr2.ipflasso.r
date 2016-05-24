### Cross-validated integrative lasso with cross-validated penalty factors
###
### Copyright 2015-07 Anne-Laure Boulesteix 
###
### Runs cvr.glmnet giving different penalty factors to predictors from different blocks
###
###
### This file is part of the `ipflasso' library for R and related languages.
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
 


cvr2.ipflasso<-function(X,Y,family,type.measure,standardize=TRUE,alpha=1,blocks,pflist,nfolds,ncv,nzeromax=+Inf,plot=FALSE)
{
M<-length(blocks)
nw<-length(pflist)

if (!setequal(M,sapply(pflist,FUN=length)))
 {
 stop("The length of the entries of argument pflist must equal the number of blocks.")
 }

ulblocks<-as.numeric(unlist(blocks))
if (!setequal(ulblocks,1:ncol(X)))
 {
 stop("Each predictor should be included in exactly one block.")
 }

if (family=="gaussian")
 {
 if (type.measure!="mse")
  warning("type.measure is set to mse.")
 type.measure<-"mse"
 }

if (family=="cox")
 {
 if (type.measure!="deviance")
  warning("type.measure is set to partial likelihood.")
 }

if (family=="binomial"&!is.element(type.measure,c("auc","class")))
 {
 warning("type.measure is set to class")
 type.measure<-"class"
 }



a<-list()
cvmin<-+Inf


for (j in 1:nw)
 {
 a[[j]]<-cvr.ipflasso(Y=Y,X=X,family=family,type.measure=type.measure,standardize=standardize,alpha=alpha,blocks=blocks,pf=pflist[[j]],nfolds=nfolds,ncv=ncv)
 allowedindices<-which(as.numeric(a[[j]]$nzero)<=nzeromax)
 if (type.measure=="auc")
  {
  ajcvm<--a[[j]]$cvm
  }
 if (type.measure!="auc")
  {
  ajcvm<-a[[j]]$cvm
  }
 mincvmj<-min(ajcvm[allowedindices]) 
 if (mincvmj<cvmin)
  {
  ind.bestpf<-j
  ind.bestlambda<-which.min(ajcvm[allowedindices])[1]
  bestlambda<-a[[j]]$lambda[ind.bestlambda]
  cvmin<-mincvmj
  }
 }

if (plot==TRUE)
 {
 par(mfrow=c(1,2))
 plot(a[[1]]$nzero,a[[1]]$cvm,type="l",ylim=c(0,1),ylab=type.measure,xlab="total number of included variables")
 if (nw>1)
 {
 for (j in 2:nw)
  {
  points(a[[j]]$nzero,a[[j]]$cvm,type="l",col=j)
  }
 abline(v=nzeromax,lty=2)
 }
 pfnames<-c()
 for (j in 1:nw)
  {
  pfnames<-c(pfnames,paste(pflist[[j]],collapse="-"))
  }
 legend(col=1:nw,lty=1,legend=pfnames,y=0.95,x=max(a[[1]]$nzero)/3) 
 
 plot(a[[ind.bestpf]]$nzero,apply(a[[ind.bestpf]]$coeff[blocks[[1]]+1,],MARGIN=2,FUN=function(x)sum(x!=0)),ylim=c(0,max(a[[ind.bestpf]]$nzero)),col=ind.bestpf,pch=2,xlab="total number of included variables",ylab="number of variables")
 
 for (b in 2:length(blocks))
  {
  points(a[[ind.bestpf]]$nzero,apply(a[[ind.bestpf]]$coeff[blocks[[b]]+1,],MARGIN=2,FUN=function(x)sum(x!=0)),pch=b+1,col=ind.bestpf)
  }
 abline(a=0,b=1)
 legend(pch=1+(1:length(blocks)),legend=paste("block",1:length(blocks)),x=0,y=max(a[[ind.bestpf]]$nzero),col=ind.bestpf)
 }

coeff<-a[[ind.bestpf]]$coeff

return(list(coeff=coeff,ind.bestlambda=ind.bestlambda,bestlambda=bestlambda,ind.bestpf=ind.bestpf,a=a,family=family))
}
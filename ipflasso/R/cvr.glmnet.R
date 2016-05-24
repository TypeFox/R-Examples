###   Repeating cv.glmnet
###
### Copyright 2015-07 Anne-Laure Boulesteix 
###
### 
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
 


cvr.glmnet<-function(X,Y,family,standardize=TRUE,alpha=1,nfolds,ncv,type.measure,...)
{
set.seed(1)
a1<-cv.glmnet(x=X,y=Y,family=family,standardize=standardize,nfolds=nfolds,type.measure=type.measure,alpha=alpha,...)
lambda1<-a1$lambda
cvm<-a1$cvm

if (ncv==1)
{
return(list(coeff=rbind(as.numeric(a1$glmnet.fit$a0),
                        as.matrix(a1$glmnet.fit$beta)),lambda=lambda1,cvm=cvm))
}

else
{
lambda<-lambda1
cvm<-matrix(cvm,nrow=1)
for (i in 2:ncv)
 {
 set.seed(i)
 a<-cv.glmnet(x=X,y=Y,family=family,standardize=standardize,nfolds=nfolds,type.measure=type.measure,lambda=lambda1,alpha=alpha,...)
 newlambda<-intersect(lambda,a$lambda)
 cvm<-rbind(cvm[,is.element(lambda,newlambda)],a$cvm[is.element(a$lambda,newlambda)])
 lambda<-newlambda
 }

coeff<-rbind(as.numeric(a1$glmnet.fit$a0)[1:length(lambda)],
             as.matrix(a1$glmnet.fit$beta)[,1:length(lambda)])

return(list(coeff=coeff,lambda=lambda,cvm=apply(cvm,MARGIN=2,FUN=mean)))
}
}
###   Cross-validated lasso with fixed penalty factors
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
 



cvr.ipflasso<-function(X,Y,family,type.measure,standardize=TRUE,alpha=1,blocks,pf,nfolds,ncv)
{
M<-length(blocks)
if (M!=length(pf))
 {
 stop("blocks and pf must have the same length.")
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



if (any(pf<=0))
 {
 stop("pf should have positive entries.")
 }


pfvector<-numeric(ncol(X))
for (m in 1:M)
 {
 pfvector[blocks[[m]]]<-pf[m]
 }


a<-cvr.glmnet(X=X,Y=Y,family=family,standardize=standardize,alpha=alpha,nfolds=nfolds,ncv=ncv,type.measure=type.measure,penalty.factor=pfvector)
coeff<-a$coeff

if (family!="cox")
 {
 rownames(coeff)[1]<-"intercept"
 }


if (type.measure=="auc")
 {
 ind.bestlambda<-which.max(a$cvm)
 }

if (type.measure!="auc")
 {
 ind.bestlambda<-which.min(a$cvm)
 }

nzero<-apply(coeff[-1,],FUN=function(x) return(sum(x!=0)),MARGIN=2) 

return(list(coeff=coeff,ind.bestlambda=ind.bestlambda,lambda=a$lambda,cvm=a$cvm,nzero=nzero,family=family))
}

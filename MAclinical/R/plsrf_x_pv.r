###    Class prediction based on PLS dimension reduction (without pre-validation) and random forests 
### 	with microarray data and clinical parameters
###
### Copyright 2007-11 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `MAclinical' library for R and related languages.
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


plsrf_x_pv<-function(Xlearn,Zlearn=NULL,Ylearn,Xtest,Ztest=NULL,ncomp=0:3,ordered=NULL,nbgene=NULL,fold=10,...)
{
ncomp<-ncomp[ncomp!=0]
n<-length(Ylearn)
ncv<-n/fold
cv<-floor((0:(n-1))/ncv)+1
Ylearn<-as.numeric(factor(Ylearn))-1
nlearn<-length(Ylearn)
p<-ncol(Xlearn)

if (is.null(ordered)&is.null(nbgene))
 {
 ordered<-1:p
 }

if (is.null(ordered)&!is.null(nbgene))
 {
 ordered<-order(abs(studentt.stat(X=Xlearn,L=Ylearn+1)),decreasing=TRUE)
 }

if (is.null(nbgene))
 {
 nbgene<-p
 }


XXlearn<-matrix(0,nlearn,max(ncomp))
for (i in 1:fold)
 { 
 print(paste("fold",i))
 if (nbgene<ncol(Xlearn))
  {
  orderedi<-order(abs(studentt.stat(X=Xlearn[cv!=i,],L=Ylearn[cv!=i]+1)),decreasing=TRUE)
  }
 else
  {
  orderedi<-1:ncol(Xlearn)
  }
 output.plsi<-pls.regression(Xlearn[cv!=i,orderedi[1:nbgene]],Ylearn[cv!=i]+1,ncomp=max(ncomp))
 XXlearn[cv==i,]<-scale(matrix(Xlearn[cv==i,orderedi[1:nbgene]],nrow=sum(cv==i)),scale=FALSE,center=output.plsi$meanX)%*%output.plsi$R[,1:max(ncomp)]
 }

output.pls<-pls.regression(Xlearn[,ordered[1:nbgene]],Ylearn+1,ncomp=max(ncomp))
XXtest<-matrix(scale(Xtest[,ordered[1:nbgene]],scale=FALSE,center=output.pls$meanX)%*%output.pls$R[,1:max(ncomp)],nrow=nrow(Xtest))
output.forest<-list()
OOB<-numeric(length(ncomp))


for (i in 1:length(ncomp))
  {
  data.learn<-data.frame(XXlearn[,1:ncomp[i]],y=factor(Ylearn))
  data.test<-data.frame(XXtest[,1:ncomp[i]])
  names(data.learn)<-c(sapply(as.list(1:ncomp[i]),FUN=paste,".comp"),"y")
  names(data.test)<-sapply(as.list(1:ncomp[i]),FUN=paste,".comp")
  output.forest[[i]]<-cforest(formula=y~.,data=data.learn,controls=cforest_control(ntree=200,mincriterion=qnorm(0.5),mtry=floor(sqrt(ncol(data.learn)-1)),replace=FALSE))
  OOB[i]<-sum(predict(output.forest[[i]],OOB=TRUE)!=Ylearn)/nlearn
  }


best<-which.min(OOB)
bestncomp<-ncomp[best]
output.forest<-output.forest[[best]]
importance<-varimp(output.forest)
  
newdata<-data.frame(data.test[,1:bestncomp])
names(newdata)<-sapply(as.list(1:bestncomp),FUN=paste,".comp")
prediction<-as.numeric(predict(object=output.forest,newdata=newdata))-1

return(list(prediction=prediction,importance=importance,bestncomp=bestncomp,OOB=OOB))
}

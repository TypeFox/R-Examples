#########################################################################
#########################################################################
#########################################################################
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License v.3 as published by
# the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# See <http://www.gnu.org/licenses/>
#    
#   
#########################################################################
#########################################################################
######################################################################### 

KODAMA = function(data,
                  M         = 100,
                  Tcycle    = 20,
                  FUN_VAR   = function(x) {ceiling(ncol(x))},
                  FUN_SAM   = function(x) {ceiling(nrow(x)*0.75)},
                  bagging   = FALSE,
                  FUN       = KNN.CV,
                  f.par     = list(kn=10),
                  W         = NULL,
                  constrain = NULL,
                  fix       = rep(FALSE,nrow(data)),
                  epsilon   = 0.05,
                  shake     = FALSE){

   n=M
   m=Tcycle
   startv=W
   data = as.matrix(data)
   nva=ncol(data)
   nsa=nrow(data)
   ma=matrix(0,ncol=nsa,nrow=nsa)
   normalization=matrix(0,ncol=nsa,nrow=nsa)   
   FUN_VAR=FUN_VAR(data)
   FUN_SAM=FUN_SAM(data)
   vect_acc=matrix(nrow=n,ncol=Tcycle)
   accu=NULL
   whF=which(!fix)
   whT=which(fix)   
   FUN_SAM=FUN_SAM-length(whT)
   pb <- txtProgressBar(min = 0, max = n, style = 1)
   for(k in 1:n){
      setTxtProgressBar(pb, k)
      sva=sample(nva,FUN_VAR, FALSE, NULL)    #.Internal(sample(nva,FUN_VAR, FALSE, NULL))      
      ssa=c(whT,sample(whF,FUN_SAM, bagging, NULL))
      x=data[ssa,sva]
      xva=ncol(x)
      xsa=nrow(x)
      if(is.vector(startv))    clbest=startv[ssa]
      if(is.function(startv))  {
         clbest=startv(x)
      }else{ 
         if(is.null(startv)) {    
            clbest=1:xsa
         }else{ 
            if(any(is.na(startv))) {
               unw=unique(startv)
               ull=length(unw)
               unw=unw[-which(is.na(unw))]     
               ghg=is.na(startv)
               startv[ghg]=(1:(sum(ghg)+length(unw)))[-unw]
               clbest=startv[ssa]
            }
         }
      }
      yatta=core(x,clbest, Tcycle=Tcycle,FUN=FUN,f.par=f.par,constrain=constrain[ssa],fix=fix[ssa],shake)
      clbest=yatta$c
      accu=yatta$a
      knk=length(yatta$v)
      vect_acc[k,1:knk]=yatta$v
      uni=unique(clbest)
      nun=length(uni)
      for(ii in 1: nun)
         ma[ssa[clbest==uni[ii]],ssa[clbest==uni[ii]]] = ma[ssa[clbest==uni[ii]],ssa[clbest==uni[ii]]]+1
      normalization[ssa,ssa]=normalization[ssa,ssa]+1
   }
   close(pb)
   ma=ma/normalization
   Edist=as.matrix(dist(data))
   ma[ma<epsilon]=0
   mam=(1/ma)*Edist
   mam=allShortestPaths(mam)$length
   prox=Edist/mam
   diag(prox)=1
   prox[is.na(prox)]=0
   mam[is.na(mam)]=max(mam,na.rm=T)
   mam=as.dist(mam)
   return(list(dissimilarity=mam,acc=accu,proximity=prox,v=vect_acc))
}

#############################################################                          

core = function(x,              #matrix
                clbest,         #starting vector
                Tcycle=20,           #number of cycles
                FUN=KNN.CV,
                f.par=list(kn=10),
                constrain=NULL,
                fix=NULL,
                shake=FALSE)
{
   m=Tcycle
   if(is.null(constrain))   constrain=1:length(clbest)
   if(is.null(fix))         fix=rep(FALSE,length(clbest))
   xsa=nrow(x)
   nconc=length(unique(constrain))   
   clbest=as.factor(clbest)  
   cvpred=clbest
   cvpredbest=do.call(FUN, c(list(x,clbest,constrain),f.par))    
   levels(cvpredbest)=levels(clbest)  
   if(shake==FALSE){
      accbest=mean(as.numeric(clbest==cvpredbest))   
   }else{
      accbest=0	
   }
   success=FALSE
   j=0
   vect_acc=NULL
   while(j<m & !success){
      j=j+1
      cl=clbest
      temp=sample(nconc, 1, FALSE, NULL)#.Internal(sample(nconc, 1, FALSE, NULL))
      ss=sample(nconc, temp, FALSE, NULL)#.Internal(sample(nconc, temp, FALSE, NULL))
      for(k in ss) {       
         sele=(constrain==k & !fix)     
         if(sum(sele!=0)) {
            soso=sample(unique(cvpredbest[sele]),1)
            cl[sele]=soso     
         } 
      }
      cvpred=do.call(FUN, c(list(x,(cl),constrain),f.par))
      levels(cvpred)=levels(cl)
      accTOT=mean(as.numeric(cl==cvpred))
      if(accTOT>accbest){
         cvpredbest=cvpred
         clbest=cl       
         accbest=accTOT
      }
      vect_acc[j]=accbest
      if(accTOT==1)   success=TRUE
   }
   return(list(c=clbest,a=accbest,v=vect_acc))
}

#############################################################                          

PCA.CA.KNN.CV = function(x,cl,constrain,kn=10,variance=0.9){
   cl = as.factor(cl)
   uni= unique(cl)
   xsa=nrow(x)
   if(length(uni)==1) return(rep(uni,xsa))
   pr=prcomp(x)
   va=pr$sdev^2
   va=va/sum(va)
   cum=which(cumsum(va)<variance)
   if(length(cum)<2) cum=1:2
   x=pr$x[,cum]
   ncomp=ncol(x)
   ca.out <- cancor(x, transformy(cl))
   nn <- length(ca.out$cor)
   mm <- ncol(ca.out$xcoef)
   if (nn == 1 && mm > 1)  nn <- nn + 1
   if (mm < ncomp)   ca.out$xcoef = rbind(ca.out$xcoef,matrix(0,ncol=mm,nrow=ncomp-mm))
   x <- x %*% ca.out$xcoef[,1:nn,drop=F]
   cvpred=cl
   fold = kfold(constrain) 
   xdist=knn.dist(x)
   for(i in 1:10){
      ff=fold==i
      w1=which(ff)
      w9=which(!ff)
      cvpred[w1]=knn.predict(train=w9,test=w1,cl,xdist,k=kn)
   }
   cvpred
}


#############################################################    

KNN.CV = function(x,cl,constrain,kn=10){
   xsa=nrow(x)
   uni = unique(cl)
   cl=as.factor(cl)
   if(length(uni)==1) return(rep(uni,xsa))
   cvpred = cl
   fold = kfold(constrain)  
   xdist=knn.dist(x)
   for(i in 1:10){
      ff=fold==i
      w1=which(ff)
      w9=which(!ff)
      cvpred[w1]=knn.predict(train=w9,test=w1,cl,xdist,k=kn)
   }
   cvpred
}

#############################################################                          
 
PLS.SVM.CV = function(x,cl,constrain,ncomp=5, ...){
   xsa=nrow(x)
   pls.svm.predict = function (train,test,cl,data,ncomp=5, ...){
      Xtrain=data[train,]
      Xtest=data[test,,drop=F]
      Ytrain=cl[train]
      uni=unique(Ytrain)
      if(length(uni)==1) return(rep(uni,length(test)))
      ntrain <- nrow(Xtrain)
      nn <- min(dim(Xtrain))-1
      if ( ncomp == 0 | ncomp > nn) ncomp <- nn
      Xtrain <- scale(Xtrain, center=T, scale=F)
      Xtest <-  scale(Xtest, center=attr(Xtrain,"scaled:center"), scale=F)
      red.out = pls.regression(Xtrain = Xtrain, transformy(Ytrain), Xtest = NULL, ncomp = ncomp)
      new.train <- matrix(red.out$T[,1:ncomp], nrow=ntrain, ncol=ncomp)
      new.test <- Xtest %*% red.out$R[,1:ncomp]

      cl.out <- svm(x=new.train, y=Ytrain,type="C-classification", ... )
      predict(object = cl.out, newdata = new.test)
   }
   cvpred=NULL
   fold = kfold(constrain)
   for(i in 1:10){
      ff=fold==i
      w1=which(ff)
      w9=which(!ff)
      if(length(unique(cl[w9]))>1){
         cvpred[w1]=pls.svm.predict(train=w9,test=w1,cl,x,ncomp=ncomp,...) 
      }else{
         cvpred[w1]=cl[w1]
      }
      
   }
   as.factor(cvpred)
} 

#########################################################

kfold = function(constrain,k=10){
   constrain=as.numeric(as.factor(constrain))
   xsa=max(constrain)
   v=sample(1:xsa)
   (v[constrain]%%k+1)
}

#############################################################

transformy = function (y) {
   y  =  as.numeric (as.factor( y ))
   n  = length(y)
   nc = max(y,na.rm=T)
   Y  = matrix(0, n, nc)
   for(k in 1:nc) 
      Y[, k] <- as.numeric(y == k)
   Y
}

#############################################################

majority <- function(x){
  x <- as.factor(x)
  n <- nlevels(x)
  votes <- rep(0,n)
  for (i in 1:length(x)) votes[as.integer(x[i])] <- votes[as.integer(x[i])]+1
  levels(x)[order(votes,decreasing=TRUE,sample(1:n,n))[1]]
  }

#########################################################

knn.probability <- function(train, test, y, dist.matrix, k=1, ties.meth="min") {

#number of predictions to make
n<-length(test)

#sort the indexes for the training and test sets
if (is.unsorted(train)) train<-sort(train)
if (is.unsorted(test)) test<-sort(test)

#only need the rows for the test data and columns
#for the training data
d<-dist.matrix[test,train]

#ensure y is a factor
y<-as.factor(y)

#only need the responses for the training data
if (length(y)>length(train)) y<-y[train]

#calculate closest neighbors and
#return aggregate response for the k closest neighbors
if (n==1) {
  d<-rank(d, ties.method = ties.meth)
  x<-classprob(y[d <= k])
  x<-data.frame(x)
  names(x)<-test
  row.names(x)<-levels(y)
  return(x)
  }
else {
  d<-t(apply(d,1,function(x) rank(x,ties.method=ties.meth)))
  x<-apply(d,1,function(x) classprob(y[x<=k]))
  row.names(x)<-levels(y)
  return(x)
  }
}

#########################################################

knn.predict <- function(train, test, y, dist.matrix, k=1,
    agg.meth=if (is.factor(y)) "majority" else "mean",
    ties.meth="min") {

#number of predictions to make
n<-length(test)

#sort the indexes for the training and test sets
if (is.unsorted(train)) train<-sort(train)
if (is.unsorted(test)) test<-sort(test)

#only need the rows for the test data and columns
#for the training data
d<-dist.matrix[test,train]

#only need the responses for the training data
if (length(y)>length(train)) y<-y[train]

#calculate closest neighbors and
#return aggregate response for the k closest neighbors
if (n==1) {
  d <- rank(d, ties.method = ties.meth)
  x <- apply(data.frame(y[d <= k]), 2, agg.meth)
  names(x) <- test
  return(x)
  }
else {
  d<-t(apply(d,1,function(x) rank(x,ties.method=ties.meth)))
  apply(d,1,function(x) apply(data.frame(y[x<=k]),2,agg.meth))
  }
}

#########################################################

knn.dist <- function(x, dist.meth="euclidean", p=2) {
   #create a distance matrix using all values in the data
   d<-as.matrix(dist(x,dist.meth,p))
   #fix for some small high persision errors
   round(d,digits=15)
}

#########################################################

classprob <- function(x){
  x <- as.factor(x)
  n <- nlevels(x)
  votes <- rep(0, n)
  for (i in 1:length(x)) votes[as.integer(x[i])] <- votes[as.integer(x[i])]+1
  votes/length(x)
}

#########################################################

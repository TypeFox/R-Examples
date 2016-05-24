align.missing <-
function(X,nlandmarks){
 
  
  
  ##########################################
  ########### internal functions############
  ##########################################  
  
  
  format.array<-function(dataset,nlandmarks){
    dataset<-as.matrix(dataset)
    nspecimen<-nrow(dataset)/nlandmarks
    start<-seq(from=1,to=nrow(dataset),by=nlandmarks)
    sparray<-array(dataset,dim<-c(nlandmarks,2,nspecimen))
    for (k in 1:nspecimen){
      x<-start[k]
      y<-x+nlandmarks-1
      single<-dataset[x:y,]
      sparray[,,k]<-single
    }
    return(sparray)
  }
  
  
  
  complete.specimens<-function(dataset,nlandmarks){
    base<-c(1,1)
    included<-base
    excluded<-base
    nspecimen<-nrow(dataset)/nlandmarks
    start<-seq(from=1,to=nrow(dataset),by=nlandmarks)
    nsp<-1:nspecimen
    for (k in 1:nspecimen){
      x<-start[k]
      y<-x+nlandmarks-1
      single<-dataset[x:y,]
      reduced<-na.omit(single)
      rows<-nrow(reduced)
      if (rows==nlandmarks){included<-rbind(included,single)}
      else {excluded<-rbind(excluded,single)}
    }
    end<-nrow(included)
    included<-included[2:end,]
    return (included)
  }
  
  ########## specimen alignment ##########
  
  char.com<-matrix(ncol=ncol(X),nrow=nrow(X))
  xvalues<-X[,1]
  nspecimen<-nrow(X)/nlandmarks
  xmat<-matrix(xvalues,ncol=nlandmarks,nrow=nspecimen,byrow=TRUE)
  missguide<-ifelse((is.na(xmat))==TRUE,1,0)
  whichpoints<-missguide*matrix(1:nlandmarks,ncol=nlandmarks,nrow=nspecimen,byrow=TRUE)
  eachmissing<-apply(missguide,1,sum)
  incompletes<-setdiff(ifelse(eachmissing==0,0,1)*1:nspecimen,0)
  
  completes<-complete.specimens(X,nlandmarks)
  GPA.com<-procGPA(format.array(completes,nlandmarks),pcaoutput=FALSE,distances=FALSE)
  aligned.com<-GPA.com$rotated
  mean.aligned<-GPA.com$mshape
  
  aligned.mat<-aligned.com[,,1]
  for(z in 2:dim(aligned.com)[3]){
    aligned.mat<-rbind(aligned.mat,aligned.com[,,z])
  }
  
  for (i in 1:length(incompletes)){
    current<-incompletes[i]  
    cur.spec<-X[((current-1)*nlandmarks+1):(current*nlandmarks),]
    cur.miss<-setdiff(whichpoints[current,],0)
    cur.spec.miss<-cur.spec[-cur.miss,]
    mean.miss<-mean.aligned[-cur.miss,]
    OPA<-procOPA(mean.miss,as.matrix(cur.spec.miss))
    transpose<-mean.miss[1,]-OPA$Ahat[1,]
    new.specX<-OPA$Bhat[,1]+transpose[1]
    new.specY<-OPA$Bhat[,2]+transpose[2]
    new.spec<-cbind(new.specX,new.specY)
    new.spec.in<-new.spec
    for (m in 1:length(cur.miss)){
      insert<-cur.miss[m]
      new.spec.in<-insertRow(new.spec.in,insert)
    }
    starts<-(current-1)*nlandmarks+1
    stops<-current*nlandmarks
    char.com[starts:stops,]<-new.spec.in
  }
  
  completed<-setdiff(1:nspecimen,incompletes)
  for (f in 1:length(completed)){
    com.n<-completed[f]
    com.spec<-aligned.com[,,f]
    starts<-(com.n-1)*nlandmarks+1
    stops<-com.n*nlandmarks
    char.com[starts:stops,]<-com.spec
  }
  
  return(char.com)
  
  
  
  
  
}

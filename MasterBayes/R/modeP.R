"modeP"<-function(postP, threshold=0, marginal=FALSE, USasNA=TRUE){

    post_prob2=NULL
  if(is.list(postP)){
    ped<-matrix(NA, length(postP), 3)
    ped[,1]<-names(postP)
    if(marginal){
      ped[,2]<-unlist(lapply(postP, function(x){colnames(x)[which.max(colSums(x))]})) 
      post_prob<-unlist(lapply(postP, function(x){colSums(x)[which.max(colSums(x))]/sum(x)})) 
      ped[,2][which(post_prob<threshold)]<-NA
      if(USasNA){
        ped[,2][which(ped[,2]=="us")]<-NA
      }
      ped[,3]<-unlist(lapply(postP, function(x){rownames(x)[which.max(rowSums(x))]})) 
      post_prob2<-unlist(lapply(postP, function(x){rowSums(x)[which.max(rowSums(x))]/sum(x)})) 
      ped[,3][which(post_prob2<threshold)]<-NA
      if(USasNA){
        ped[,3][which(ped[,3]=="us")]<-NA
      }
    }else{
      ped[,2]<-unlist(lapply(postP, function(x){colnames(x)[which(x==max(x), arr.ind=TRUE)[1,][2]]})) 
      ped[,3]<-unlist(lapply(postP, function(x){rownames(x)[which(x==max(x), arr.ind=TRUE)[1,][1]]})) 
      post_prob<-unlist(lapply(postP, function(x){x[which.max(x)]/sum(x)})) 
      ped[,2][which(post_prob<threshold)]<-NA
      ped[,3][which(post_prob<threshold)]<-NA
      if(USasNA){
        ped[,2][which(ped[,2]=="us")]<-NA
        ped[,3][which(ped[,3]=="us")]<-NA
      }
    }
  }else{
    ped<-matrix(NA, dim(postP)[1], 3)
    lpost<-dim(postP)[2]
    ped[,1]<-rownames(postP)
    if(marginal){
       postP1<-apply(postP, 1, function(x){table(x[seq(1,lpost,2)])/(lpost/2)})
       if(is.list(postP1)==FALSE){
         ped[,2]<-apply(postP, 1, function(x){x[1]})
         post_prob<-rep(1, dim(ped)[1])
       }else{
         ped[,2]<-unlist(lapply(postP1, function(x){names(x)[which.max(x)]}))
         post_prob<-unlist(lapply(postP1, function(x){max(x)}))
       }
       ped[,2][which(post_prob<threshold)]<-NA
       if(USasNA){
         ped[,2][which(ped[,2]=="us")]<-NA
       }       
       postP2<-apply(postP, 1, function(x){table(x[seq(2,lpost,2)])/(lpost/2)})
       if(is.list(postP2)==FALSE){
         ped[,3]<-apply(postP, 1, function(x){x[2]})
         post_prob2<-rep(1, dim(ped)[1])
       }else{
         ped[,3]<-unlist(lapply(postP2, function(x){names(x)[which.max(x)]}))
         post_prob2<-unlist(lapply(postP2, function(x){max(x)}))
       }
       ped[,3][which(post_prob2<threshold)]<-NA
       if(USasNA){
         ped[,3][which(ped[,3]=="us")]<-NA
       }
    }else{
      postP<-apply(postP, 1, function(x){table(paste(x[seq(1,lpost,2)], x[seq(2,lpost,2)]))})
      ped[,2]<-unlist(lapply(postP, function(x){strsplit(names(x)[which.max(x)], " ")[[1]][1]})) 
      ped[,3]<-unlist(lapply(postP, function(x){strsplit(names(x)[which.max(x)], " ")[[1]][2]}))
      post_prob<-unlist(lapply(postP, function(x){x[which.max(x)]/sum(x)})) 
      ped[,2:3][which(post_prob<threshold),]<-NA
      if(USasNA){
        ped[,2][which(ped[,2]=="us")]<-NA
        ped[,3][which(ped[,3]=="us")]<-NA
      }
    }
  }
  list(P=ped, prob=as.vector(post_prob), prob.male=as.vector(post_prob2))
}
 

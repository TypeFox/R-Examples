"post.pairs"<-function(postP, threshold=0, rel="PO"){

  if(is.list(postP)){
    stop("P must be samples from the JOINT posterior of parentage")
  }else{
    post<-matrix(0, dim(postP)[1], dim(postP)[1])
    if(rel=="PO"){
      for(i in 1:(dim(postP)[2]/2)){
        post<-post+outer(rownames(postP), postP[,i*2-1], "==")
        post<-post+outer(rownames(postP), postP[,i*2], "==")
      }
      post<-post+t(post)
    }else{
      dummyD<-paste("D", letters, 1:dim(postP)[1], letters)
      dummyS<-paste("S", letters, 1:dim(postP)[1], letters)
      for(i in 1:(dim(postP)[2]/2)){
        usv<-which(postP[,i*2-1]=="us")
        postP[,i*2-1][usv]<-dummyD[1:length(usv)]
        usv<-which(postP[,i*2]=="us")
        postP[,i*2][usv]<-dummyS[1:length(usv)]
        op<-outer(postP[,i*2-1], postP[,i*2-1], "==")+outer(postP[,i*2], postP[,i*2], "==")
        diag(op)<-0
        if(rel=="S"){
          op<-(op>0)
        }
        if(rel=="FS"){
          op<-(op==2)
        }
        if(rel=="HS"){
          op<-(op==1)
        }
        post<-post+op
      }
    }
  }
  rownames(post)<-rownames(postP)
  post[which(upper.tri(post)==TRUE)]<-0
  P1<-rownames(post)[which(post>(dim(postP)[2]*threshold/2), arr.ind = TRUE)[,1]]
  P2<-rownames(post)[which(post>(dim(postP)[2]*threshold/2), arr.ind = TRUE)[,2]]
  sig.pairs<-cbind(P1, P2)
  list(P=sig.pairs, prob=(2*post/dim(postP)[2])[which(post>(dim(postP)[2]*threshold/2))])
}



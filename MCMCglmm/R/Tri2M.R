 
"Tri2M"<-function(x,lower.tri=TRUE, reverse=TRUE, diag=TRUE){

  if(is.null(dim(x))==FALSE){
    if(dim(x)[1]==dim(x)[2]){
      if(dim(x)[1]==1){
        return(x)
      }else{
        if(lower.tri){
          return(x[which(lower.tri(x, diag=diag)==TRUE)])
        }else{
          return(x[which(upper.tri(x, diag=diag)==TRUE)])
        }
      }
    }else{
      dimM<-sqrt(length(x)*2+0.25)-diag*0.5+(!diag)*0.5
 
      M<-matrix(0,dimM,dimM)
      if(is.null(names(x))==FALSE){
        colnames(M)<-diag(Tri2M(names(x)))
        rownames(M)<-diag(Tri2M(names(x)))
      }
      if(lower.tri){
        M[which(lower.tri(M, diag=diag)==TRUE)]<-x
        if(reverse==TRUE){   
          M[which(upper.tri(M)==TRUE)]<-t(M)[which(upper.tri(M)==TRUE)]
        }
      }else{
        M[which(upper.tri(M, diag=diag)==TRUE)]<-x
        if(reverse==TRUE){
          M[which(lower.tri(M)==TRUE)]<-t(M)[which(lower.tri(M)==TRUE)]
        }
      }
      return(M)
    }
  }else{
      dimM<-sqrt(length(x)*2+0.25)-diag*0.5+(!diag)*0.5
 
      M<-matrix(0,dimM,dimM)
      if(is.null(names(x))==FALSE){
        colnames(M)<-diag(Tri2M(names(x)))
        rownames(M)<-diag(Tri2M(names(x)))
      }
      if(lower.tri){
        M[which(lower.tri(M, diag=diag)==TRUE)]<-x
        if(reverse==TRUE){   
          M[which(upper.tri(M)==TRUE)]<-t(M)[which(upper.tri(M)==TRUE)]
        }
      }else{
        M[which(upper.tri(M, diag=diag)==TRUE)]<-x
        if(reverse==TRUE){
          M[which(lower.tri(M)==TRUE)]<-t(M)[which(lower.tri(M)==TRUE)]
        }
      }
      return(M)
   }
}

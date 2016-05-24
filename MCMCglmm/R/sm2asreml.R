"sm2asreml"<-function(A=NULL, rownames=NULL){
  ginv<-data.frame(Row=A@i+1, Column=rep(1:length(A@p[-1]), diff(A@p)), Ainverse=A@x)
  ginv<-ginv[order(ginv$Row),]
  ginv<-ginv[which(ginv$Row>=ginv$Column),]
  if(is.null(rownames)){
    if(is.null(rownames(A))){
       stop("A must have rownames, or rownames must be provided in the call to sm2asreml")
    }else{
      attr(ginv, "rowNames")<-rownames(A)
    }
  }else{
    attr(ginv, "rowNames")<-rownames
  }
  return(ginv)
}


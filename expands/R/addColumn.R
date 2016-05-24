.addColumn<-function(M,newCol,initVal){
  if (!any(colnames(M)==newCol)){
    if(!is.null(dim(M))){
      M=matrix(cbind(M,matrix(initVal,nrow(M),1)),nrow=nrow(M),ncol=ncol(M)+1,
               dimnames = list(rownames(M), c(colnames(M),newCol)));
    }else{
      cols=names(M);
      M=c(M,initVal);	
      names(M)=c(cols,newCol);
    }
  }
  return(M);
}

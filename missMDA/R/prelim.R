prelim<-function(res.mi,X){
  if(any(c("MIMCA","MIPCA")%in%class(res.mi))){
  longformat<-rbind(X,do.call(rbind,res.mi$res.MI))
  longformat<-cbind(.imp=rep(0:length(res.mi$res.MI),each=nrow(X)),
                    .id=rep(1:nrow(X),(length(res.mi$res.MI)+1)),
                    longformat)
  rownames(longformat)<-NULL
  imp.mids<-as.mids(longformat)
  }else{
    stop("prelim requires as input an object of class MIPCA or MIMCA.")
    }
  return(imp.mids)
}

"one2two" <-
function(M,clu=NULL){
  if(!is.null(clu)){
    if(mode(clu)=="list"){
      n<-sapply(clu,FUN=length)
      newM<-M[1:n[1],(n[1]+1):sum(n[1:2])]
    } else stop("For now clu must be supplied in form of a list (one component for each mode)")
  } else stop("For now clu must be supplied in form of a list (one component for each mode)")
  return(list(M=newM,clu=clu))
}


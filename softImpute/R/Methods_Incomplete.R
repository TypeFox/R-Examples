setClass("Incomplete","dgCMatrix")
as.matrix.Incomplete=function(x,...){
  x=as(x,"dgTMatrix")
  i=x@i
  j=x@j
  out=as.matrix(x)
  out[]=NA
  nrow=dim(x)[1]
  out[i+1 +nrow*j]=x@x
  out
}
setMethod("as.matrix","Incomplete",as.matrix.Incomplete)
complete.matrix=function(x,object,unscale=TRUE){
  nas=is.na(x)
  if(!any(nas))return(x)
  i=row(x)[nas]
  j=col(x)[nas]
  x[nas]=impute(object,i,j,unscale=unscale)
  x
}
scomplete = function(x,object,unscale=TRUE){
###x is class Incomplete
    x=as.matrix(x)##Fills in NAs where the zeros are
    complete(x,object,unscale=unscale)
  }

setGeneric("complete",complete.matrix)
setMethod("complete","Incomplete",scomplete)

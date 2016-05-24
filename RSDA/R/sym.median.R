sym.median <-
function(sym.var,method=c('centers','interval','histogram'),na.rm=FALSE, ...) {
  method<-match.arg(method)
  if(method=='centers') {
    if(sym.var$var.type=='$C')
      return(median(sym.var$var.data.vector,na.rm))  
    if(sym.var$var.type=='$I') 
      return(median(sym.var$var.data.vector[,1]+sym.var$var.data.vector[,2])/2)
    else 
      stop("Impossible to compute the median for this type of variable with this method")    
  }
  if(method=='interval') {
    if(sym.var$var.type=='$I') 
      return(sapply(sym.var$var.data.vector,median))
    else 
      stop("Impossible to compute the median for this type of variable with this method")    
  }  
  if(method=='histogram') {
    if(sym.var$var.type=='$H')
      return(sapply(sym.var$var.data.vector,median))
    else 
      stop("Impossible to compute the median for this type of variable with this method")    
  }
}

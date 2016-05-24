sym.mean <-
function(sym.var,method=c('centers','interval','histogram'),trim=0,na.rm=FALSE, ...) {
  method<-match.arg(method)
  if(method=='centers') {
    if(sym.var$var.type=='$C')
      return(mean(sym.var$var.data.vector,trim,na.rm))  
    if(sym.var$var.type=='$I') 
      return(mean((sym.var$var.data.vector[,1]+sym.var$var.data.vector[,2])/2))
    else 
      stop("Impossible to compute the mean for this type of variable with this method")    
  }
  if(method=='interval') {
    if(sym.var$var.type=='$I') 
      return(colMeans(sym.var$var.data.vector))
    else 
      stop("Impossible to compute the mean for this type of variable with this method")    
  }  
  if(method=='histogram') {
    if(sym.var$var.type=='$H')
      return(colMeans(sym.var$var.data.vector))
    else 
      stop("Impossible to compute the mean for this type of variable with this method")    
  }  
}

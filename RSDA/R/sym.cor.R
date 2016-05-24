sym.cor <-
function(sym.var.x,sym.var.y,method=c('centers','interval','billard','histogram'),na.rm=FALSE, ...) {
  method<-match.arg(method)
  if(method=='centers') {
    if((sym.var.x$var.type=='$C')&&(sym.var.y$var.type=='$C'))
      return(cor(sym.var.x$var.data.vector,sym.var.y$var.data.vector))  
    if((sym.var.x$var.type=='$I')&&(sym.var.y$var.type=='$I'))
      return(cor((sym.var.x$var.data.vector[,1]+sym.var.x$var.data.vector[,2])/2,
                 (sym.var.y$var.data.vector[,1]+sym.var.y$var.data.vector[,2])/2))
    else 
      stop("Impossible to compute the Standard Deviation for this type of variable with this method")    
  }
  if(method=='billard') {
    if((sym.var.x$var.type=='$I')&&(sym.var.y$var.type=='$I'))
      return(sym.cov(sym.var.x,sym.var.y,method='billard')/(sym.sd(sym.var.x,method='billard')*sym.sd(sym.var.y,method='billard')))
  }
}

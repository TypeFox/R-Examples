sym.lm.bi <-
function(sym.var.x,sym.var.y,method=c('mid-points','tops',
                                              'inf-sup','billard')) {
  method<-match.arg(method)
  if(((sym.var.x$var.type!='$C')||(sym.var.y$var.type!='$C'))&&
       ((sym.var.x$var.type!='$I')||(sym.var.y$var.type!='$I')))
    stop("Impossible to use lm this type of variable")    
  if(method=='mid-points') {
    if((sym.var.x$var.type=='$C')&&(sym.var.y$var.type=='$C')) {
      lm1<-lm(sym.var.y$var.data.vector~sym.var.x$var.data.vector)
    }
    if((sym.var.x$var.type=='$I')&&(sym.var.y$var.type=='$I')) {
      vx<-(sym.var.x$var.data.vector[,1]+sym.var.x$var.data.vector[,2])/2
      vy<-(sym.var.y$var.data.vector[,1]+sym.var.y$var.data.vector[,2])/2
      lm1<-lm(vy~vx)
    }
    return(lm1)
  }
  if(method=='inf-sup') {
    if((sym.var.x$var.type=='$I')&&(sym.var.y$var.type=='$I')) {
      vx<-sym.var.x$var.data.vector[,2]
      vy<-sym.var.y$var.data.vector[,1]
      lm1<-lm(vy~vx)
      vx<-sym.var.x$var.data.vector[,1]
      vy<-sym.var.y$var.data.vector[,2]
      lm2<-lm(vy~vx)
      
    }
    return(list(inf=lm1,sup=lm2))
  }  
  if(method=='tops') {
    if((sym.var.x$var.type=='$I')&&(sym.var.y$var.type=='$I')) {
      vx<-c(sym.var.x$var.data.vector[,1],sym.var.x$var.data.vector[,1],
            sym.var.x$var.data.vector[,2],sym.var.x$var.data.vector[,2])
      vy<-c(sym.var.y$var.data.vector[,1],sym.var.y$var.data.vector[,2],
            sym.var.y$var.data.vector[,1],sym.var.y$var.data.vector[,2])
      lm1<-lm(vy~vx)  
    }
    return(lm1)
  }    
  if(method=='billard') {
    if((sym.var.x$var.type=='$I')&&(sym.var.y$var.type=='$I')) {
      vx<-sym.var.x
      vy<-sym.var.y
      beta1<-sym.cov(vx,vy,method='billard')/sym.variance(vx,method='billard')
      intercept<-sym.mean(vy)-beta1*sym.mean(vx)
    }
    return(list(Intercept=intercept,Beta1=beta1))   
  }
}

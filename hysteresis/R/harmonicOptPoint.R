harmonicOptPoint <-
  function(x,y,ti){
    n <- length(x)
    Ta.lm<-lm.fit(cbind(rep(1,n),cos(ti)),x)            
    b.x<-coef(Ta.lm)[[2]]
    cx<-coef(Ta.lm)[[1]]
    Tb.lm<-lm.fit(cbind(rep(1,n),sin(ti),cos(ti)),y)
    
    b.y<-coef(Tb.lm)[[3]]
    retention<- coef(Tb.lm)[[2]]
    cy<-coef(Tb.lm)[[1]]

    var.x <- crossprod(Ta.lm$residuals)/(n-2)
    var.y <- crossprod(Tb.lm$residuals)/(n-3)
    
    times<-numeric(n)
    for (i in 1:length(x)) {
      x0<-x[i]
      y0<-y[i]
      zmin1<-optimize(harmonicspot,c(0,pi),"x0"=x0,"y0"=y0,"cx"=cx,"cy"=cy,"b.x"=b.x,"b.y"=b.y,"retention"=retention,"var.x"=var.x,"var.y"=var.y)
      zmin2<-optimize(harmonicspot,c(pi,2*pi),"x0"=x0,"y0"=y0,"cx"=cx,"cy"=cy,"b.x"=b.x,"b.y"=b.y,"retention"=retention,"var.x"=var.x,"var.y"=var.y)
      times[i]<-ifelse(zmin1$objective < zmin2$objective, zmin1, zmin2)[[1]]
    }
    pred.x<-cx +b.x*cos(times)
    pred.y<-cy +b.y*cos(times)+retention*sin(times)
    
    SSE <- crossprod(x-pred.x)+crossprod(y-pred.y)
 
    list("values"=c("b.x"=b.x,"b.y"=b.y, "cx"=cx,"cy"=cy,"retention"=retention),
           "ti"=times,"pred.x"=pred.x,"pred.y"=pred.y,"SSE"=SSE)
  }

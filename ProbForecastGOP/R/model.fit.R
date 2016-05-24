"model.fit" <-
function(par.model,emp.variog,max.dist.fit,init.val,init.var,fix.nugget){
  if(length(max.dist.fit)==0){
    l.emp.variog <- length(emp.variog$empir.variog)
    h.bin.size <- (emp.variog$bin.midpoints[l.emp.variog]-emp.variog$bin.midpoints[(l.emp.variog-1)])/2
    max.dist.fit <- (emp.variog$bin.midpoints[l.emp.variog]+h.bin.size)/(2*sqrt(2))
  }
  d <- emp.variog$bin.midpoints[0 < emp.variog$bin.midpoints & emp.variog$bin.midpoints <= max.dist.fit]
  variog <- emp.variog$empir.variog[0 < emp.variog$bin.midpoints & emp.variog$bin.midpoints <= max.dist.fit]
  w <- emp.variog$number.pairs[0 < emp.variog$bin.midpoints & emp.variog$bin.midpoints <= max.dist.fit]
  if(length(init.val)==0){
   if(length(fix.nugget)==2){
    nugget.init.val <- as.numeric(fix.nugget[2])
   }
   if(length(fix.nugget)<2){
     nugget.init.val <- lowess(emp.variog$bin.midpoints,emp.variog$empir.variog)$y[1]
   }
   var.init.val <- init.var - nugget.init.val
  }
# this is for the case the parametric model is exponential
  if(par.model=="exponential"){
    if(length(init.val)==0){
      val.r <- nugget.init.val + var.init.val*(1-exp(-1))
      r.init.val <- d[abs(variog-val.r)==min(abs(variog-val.r))]
    }
  }
# this is for the case the parametric model is spherical
  if(par.model=="spherical"){
    if(length(init.val)==0){
      r.init.val <- d[abs(variog-init.var)==min(abs(variog-init.var))]
    }
  }
# this is for the case the parametric model is gauss
  if(par.model=="gauss"){
    if(length(init.val)==0){
      val.r <- nugget.init.val + var.init.val*(1-exp(-1))
      r.init.val <- d[abs(variog-val.r)==min(abs(variog-val.r))]
    }
  }
# this is for the case the parametric model is gencauchy
  if(par.model=="gencauchy"){
    if(length(init.val)==0){
      r.init.val <- d[abs(variog-init.var)==min(abs(variog-init.var))]
      a <- 1
      b <- 2
    }
  }
# this is for the case the parametric model is whittlematern
  if(par.model=="whittlematern" | par.model=="matern"){
    if(length(init.val)==0){
      a <- 0.5
      val.r <- init.var + var.init.val*(((2^0.5)/gamma(0.5))*besselK(1,0.5))
      r.init.val <- d[abs(variog-val.r)==min(abs(variog-val.r))]
    }
  }
if(length(init.val)!=0){         
  nugget.init.val <- init.val[1]
  var.init.val <- init.val[2]
  r.init.val <- init.val[3]
  if(par.model=="whittlematern" | par.model=="matern"){
    a <- init.val[4]
  }
  if(par.model=="gencauchy"){
    a <- init.val[4]
    b <- init.val[5]
  }
}
# after determining the initial values, we can proceed to the estimation of the parameters.
if(par.model=="spherical"){
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
     f<-spher.variog.fn
     init.value <- c(var.init.val,r.init.val)
     var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0),upper=c(Inf,Inf))$par}
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
     f <- spher.variog
     init.value <- c(nugget.init.val,var.init.val,r.init.val)
     var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0,0),upper=c(Inf,Inf,Inf))$par}
  
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
     variance.var <- var.par[1]
     r.var <- var.par[2]
     par.est <- c(nugget.init.val,variance.var,r.var)
   }
        
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
     nugget.var <- var.par[1]
     variance.var <- var.par[2]
     r.var <- var.par[3]
     par.est <- c(nugget.var,variance.var,r.var)}
}

if(par.model=="whittlematern" | par.model=="matern"){
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
     f <- matern.variog.fn
     init.value <- c(var.init.val,r.init.val,a)
     var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0,0),upper=c(Inf,Inf,Inf))$par}
  
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
     f<-matern.variog
     init.value <- c(nugget.init.val,var.init.val,r.init.val,a)
     var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf))$par}
   
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
     variance.var <- var.par[1]
     r.var <- var.par[2]
     a.var <- var.par[3]
     if(a.var <= 0){
       cat(paste("Warning- Estimate of a out of parameter space:",a.var),fill=TRUE,sep="")}
     par.est <- c(nugget.init.val,variance.var,r.var,a.var)}
   
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
       nugget.var <- var.par[1]
       variance.var <- var.par[2]
       r.var <- var.par[3]
       a.var <- var.par[4]
       if(a.var <= 0){
         cat(paste("Warning- Estimate of a out of parameter space:",a.var),fill=TRUE,sep="")}
       par.est <- c(nugget.var,variance.var,r.var,a.var)}
}
if(par.model=="gencauchy"){
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
      f <- gencauchy.variog.fn
      init.value <- c(var.init.val,r.init.val,a,b)
      var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0,0,0),upper=c(Inf,Inf,2,Inf))$par}
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
      f<- gencauchy.variog
      init.value <- c(nugget.init.val,var.init.val,r.init.val,a,b)
      var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0,0,0,0),upper=c(Inf,Inf,Inf,2,Inf))$par}
       
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
      variance.var <- var.par[1]
      r.var <- var.par[2]
      a.var <- var.par[3]
      b.var <- var.par[4]
      if(a.var <= 0 | a.var > 2){
       cat(paste("Warning- Estimate of a out of parameter space:",a.var),fill=TRUE,sep="")}
      par.est <- c(nugget.init.val,variance.var,r.var,a.var,b.var)}
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
      nugget.var <- var.par[1]
      variance.var <- var.par[2]
      r.var <- var.par[3]
      a.var <- var.par[4]
      b.var <- var.par[5]
      if(a.var <= 0 | a.var > 2){
       cat(paste("Warning- Estimate of a out of parameter space:",a.var),fill=TRUE,sep="")}
       par.est <- c(nugget.var,variance.var,r.var,a.var,b.var)}
}
if(par.model=="exponential"){
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
      f<- expvariog.fn
      init.value <- c(var.init.val,r.init.val)
      var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0),upper=c(Inf,Inf))$par}
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
      f<- expvariog
      init.value <- c(nugget.init.val,var.init.val,r.init.val)
      var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0,0),upper=c(Inf,Inf,Inf))$par}
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
      variance.var <- var.par[1]
      r.var <- var.par[2]
      par.est <- c(nugget.init.val,variance.var,r.var)}
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
      nugget.var <- var.par[1]
      variance.var <- var.par[2]
      r.var <- var.par[3]
      par.est <- c(nugget.var,variance.var,r.var)}
}
if(par.model=="gauss"){
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
      f<-gauss.variog.fn
      init.value <- c(var.init.val,r.init.val)
      var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0),upper=c(Inf,Inf))$par}
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
      f <- gauss.variog
      init.value <- c(nugget.init.val,var.init.val,r.init.val)
      var.par <- optim(init.value,variog=variog,d=d,w=w,fn=f,gr=NULL,method="L-BFGS-B",lower=c(0,0,0),upper=c(Inf,Inf,Inf))$par}
   if(fix.nugget[1]=="TRUE" | fix.nugget[1]==TRUE){
      variance.var <- var.par[1]
      r.var <- var.par[2]
      par.est <- c(nugget.init.val,variance.var,r.var)}
   if(fix.nugget[1]=="FALSE" | fix.nugget[1]==FALSE){
      nugget.var <- var.par[1]
      variance.var <- var.par[2]
      r.var <- var.par[3]
      par.est <- c(nugget.var,variance.var,r.var)}
}
return(par.est)
}

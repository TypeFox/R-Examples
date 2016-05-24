fun.gld.slope.fixed.int.vary <-
function(q,fit,fit.simu,fun,param,maxit=20000,method="Brent"){

# Use a default lower and upper limit

lubound<-quantile(subset(fit.simu,select=c("(Intercept)")),c(0,1))

lbound<-min(lubound)
ubound<-max(lubound)

# First find empirical solution
x<-fit$x
y<-fit$y
k1<-apply(fit.simu,2,function(x,q) quantile(x,q),q)[1]
fit1<-fit[[3]][-c((length(fit[[3]])-3):length(fit[[3]]))]
r1<-optim(k1,function(k1,x,y,fit1,q){resid<-y-data.matrix(x)%*%c(k1,fit1[-1])
return((sum(resid<=0)/length(resid)-q)^2)
},x=x,y=y,fit1=fit1,q=q,control=list(maxit=maxit),method=method,lower=lbound,upper=ubound)

# Then find parametric solution
k2<-r1$par
r2<-optim(k2,function(k2,x,y,fit1,q){gld.fit<-fun(y-data.matrix(x)%*%c(k2,fit1[-1]))
return((pgl(0,gld.fit,param=param)-q)^2)
},x=x,y=y,fit1=fit1,q=q,control=list(maxit=maxit),method=method,lower=lbound,upper=ubound)

r.val<- setNames(c(r2$par, fit1[-1], r2$value, r2$convergence),c(names(fit1),"Objective Value","Convergence"))

return(list(r2,r.val))
}

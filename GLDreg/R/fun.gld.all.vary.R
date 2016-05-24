fun.gld.all.vary <-
function(q,fit,fit.simu,fun,param,maxit=20000,method="Nelder-Mead"){

# First find empirical solution
x<-fit$x
y<-fit$y
k1<-apply(fit.simu,2,function(x,q) quantile(x,q),q)
r1<-optim(k1,function(k1,x,y,q){resid<-y-x%*%c(k1)
return((sum(resid<=0)/length(resid)-q)^2)
},x=x,y=y,q=q,control=list(maxit=maxit),method=method)

# Then find parametric solution
k2<-r1$par
r2<-optim(k2,function(k2,x,y,q){gld.fit<-fun(y-data.matrix(x)%*%c(k2))
return((pgl(0,gld.fit,param=param)-q)^2)
},x=x,y=y,q=q,control=list(maxit=maxit),method=method)

r.val<-setNames(c(r2$par,r2$value,r2$convergence),
c(names(fit[[3]][1:(length(fit[[3]])-4)]),"Objective Value","Convergence"))

return(list(r2,r.val))

}

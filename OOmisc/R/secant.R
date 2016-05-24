secant <-
function(fun,x0,x1,eps,maxit=20,silent=FALSE){
iter=0
repeat{
iter=iter+1
f0=fun(x0)
f1=fun(x1)
x2=(x0*f1-x1*f0)/(f1-f0)
if(abs(x2-x1)<eps) break
if (iter>=maxit) stop("Maximum iteration exceeded")
x0=x1
x1=x2
if (silent==FALSE) cat("Iterasyon:",iter,", Result=",x2,fill=T)
}
if (silent==FALSE) cat("Solution:",x2,fill=T)
x2
}


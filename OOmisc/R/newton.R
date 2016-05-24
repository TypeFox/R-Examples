newton <-
function(fun,derf,x0,eps,maxit=20,silent=TRUE,tun=1){
iter=0
repeat{
iter=iter+1
x2=x0-(fun(x0)/derf(x0))/tun
if (abs(x0-x2)<=eps) break
if (iter>=maxit) stop("Maximum iteration exceeded")
x0=x2
if (silent==FALSE) cat("Iteration:",iter,", Result=",x2,fill=T)
}
if (silent==FALSE) cat("Solution:",x2,fill=T)
x2
}


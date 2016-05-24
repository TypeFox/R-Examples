damped.newton <-
function(fun,derf,x0,eps,maxit=20,damp=seq(0,40),silent=TRUE){
iter=0
repeat{
iter=iter+1
h<--fun(x0)/derf(x0)

exp1<-NULL
exp2<-fun(x0)
for (j in damp){
xnew<-x0+h/(2**j)
exp1<-c(exp1,fun(xnew))
}
ch<-exp1<exp2
ch2<-which(ch,arr.ind=TRUE)[1]-1
if (is.na(ch2)==TRUE) stop("Damping is not successful")
x2<-x0+h/(2**ch2)

if(abs(x0-x2)<eps) break
if (iter>=maxit) stop("Maximum iteration exceeded")
x0=x2
if (silent==FALSE) cat("Iteration:",iter,"; Result=",x2,fill=T)
}
if (silent==FALSE) cat("Solution:",x2,fill=T)
x2
}


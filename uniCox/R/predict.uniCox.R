
predict.uniCox=function(object,x,...){

# x is n by p

mx=object$mx
vx=object$vx
s0=object$s0
beta=object$beta
x=scale(x,center=mx,scale=F)
xs=scale(x,center=F,scale=sqrt(vx)+s0)
yhat=x%*%beta
return(yhat)
}




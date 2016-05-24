qrfit.lasso=function(X,y,tau,lambda,beta,method,maxit,toler,rho)

{
if(missing(method)){
cat("please input the method\n")
method=1
}
if(method=="mm")
return(QR.lasso.mm (X,y,tau,lambda,beta,maxit,toler))

if(method=="cd")
return(QR.lasso.cd (X,y,tau,lambda,beta,maxit,toler))

if(method=="ip")
return(QR.lasso.ip(X,y,tau,lambda))

if(method=="admm")
return(QR.lasso.admm(X,y,tau,lambda,rho,beta,maxit))

}

cqr.fit.lasso=function(X,y,tau,lambda,beta,method,maxit,toler,rho)

{
if(missing(method)){
cat("please input the method\n")
method=1
}
if(method=="mm")
return(cqr.lasso.mm(X,y,tau,lambda,beta,maxit,toler))

if(method=="cd")
return(cqr.lasso.cd(X,y,tau,lambda,beta,maxit,toler))

if(method=="admm")
return(cqr.lasso.admm(X,y,tau,lambda,rho,beta,maxit))

if(method=="ip")
cat("Sorry, the function is developing\n")
}

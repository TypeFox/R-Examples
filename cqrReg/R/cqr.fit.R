cqr.fit=function(X,y,tau,beta,method,maxit,toler,rho)

{
if(missing(method)){
cat("please input the method\n")
method=1
}
if(method=="mm")
return(cqr.mm(X,y,tau,beta,maxit,toler))

if(method=="cd")
return(cqr.cd(X,y,tau,beta,maxit,toler))

if(method=="ip")
return(cqr.ip(X,y,tau))

if(method=="admm")
return(cqr.admm(X,y,tau,rho,beta,maxit,toler))

}

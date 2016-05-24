qrfit=function(X,y,tau,beta,method,maxit,toler,rho)

{
if(missing(method)){
cat("please input the method\n")
method=1
}
if(method=="mm")
return(QR.mm(X,y,tau,beta,maxit,toler))

if(method=="cd")
return(QR.cd(X,y,tau,beta,maxit,toler))

if(method=="ip")
return(QR.ip(X,y,tau))

if(method=="admm")
return(QR.admm(X,y,tau,rho,beta,maxit))

}

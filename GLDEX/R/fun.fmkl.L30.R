"fun.fmkl.L30"<-function(k,L4){
j<-0:k

result<-integrate(function(x,k,L4) (log(x)-((1-x)^L4-1)/L4)^k,0,1,
abs.tol=1e-100,k=k,L4=L4,stop.on.error=FALSE)

if(result$message!="OK"){return(NA)}

else return(result$value)}



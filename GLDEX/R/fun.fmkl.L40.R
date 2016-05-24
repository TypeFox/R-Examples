"fun.fmkl.L40"<-function(k,L3){
j<-0:k

result<-integrate(function(x,k,L3) ((x^L3-1)/L3-log(1-x))^k,0,1,
abs.tol=1e-100,k=k,L3=L3,stop.on.error=FALSE)

if(result$message!="OK"){return(NA)}

else return(result$value)}




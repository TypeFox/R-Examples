Quantitative=function(z,p1,p2,p3,mu,sigma){
  
  if(!is.vector(z)){stop("z must be a vector.")}
  if(any(is.na(z))){stop("There are missing values in z.")}                          
  
  if((p1<0)|(p1>1)){stop("There are invalid values in p1.")}
  
  if((p2<0)|(p2>1)){stop("There are invalid values in p2.")}
  
  if((p3<0)|(p3>1)){stop("There are invalid values in p3.")}
  
  if(!is.vector(mu)){stop("mu must be a vector.")}
  if(any(is.na(mu))){stop("There are missing values in mu.")}  
  if(length(mu)!=3){stop("The length of the vector mu must be three")}
   
  if(!is.vector(sigma)){stop("sigma must be a vector.")}
  if(any(is.na(sigma))){stop("There are missing values in sigma.")}
  if(length(mu)!=length(sigma)){stop("The lengths of mu and sigma are different.")}
  if(any(sigma<0)){stop("There are invalid values in sigma.")}    
 
  r=(z-p2*mu[2]-p3*mu[3])/(p1+p2*mu[1])
  
  A=p1*(1-p1)+sigma[1]^2*p2+mu[1]^2*p2-mu[1]^2*p2^2-2*p1*p2*mu[1]
  B=2*p2*mu[1]*mu[2]-2*mu[1]*mu[2]*p2^2-2*p1*p2*mu[2]-2*mu[3]*p1*p3-2*mu[1]*mu[3]*p2*p3
  C=(sigma[2]^2+mu[2]^2)*p2+(sigma[3]^2+mu[3]^2)*p3-(mu[2]*p2+mu[3]*p3)^2
  vr=(1/(p1+p2*mu[1])^2)*(r^2*A+r*B+C)
  
  if(any(vr<0)){warning("The transformed variance estimation contains negative values.")}  
 
  out=list(TransformedVariable=r,TransformedVariance=vr)
  return(out)
}
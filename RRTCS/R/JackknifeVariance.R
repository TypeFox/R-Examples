JackknifeVariance=function(n,r,pi,strata=FALSE,cluster=FALSE,clu=NULL){
  
  if(length(n)!=1){stop("n must be a scalar.")}
  if(n<0){stop("n must be a positive number.")}
  
  if(!is.vector(r)){stop("r must be a vector.")}
  if(any(is.na(r))){stop("There are missing values in r.")}
  
  if(!is.vector(pi)){stop("pi must be a vector.")}
  if(any(is.na(pi))){stop("There are missing values in pi.")}
  if(any((pi<=0)|(pi>1))){stop("There are invalid values in pi.")}         
  if(length(pi)!=length(r)){stop("The lengths of pi and r are different.")}                 
  
  if(((strata==FALSE)&(cluster==FALSE))|((strata==TRUE)&(cluster==FALSE))){
    Je=vector()
    for(i in 1:n){
      wi=1/(pi[-i]*((n-1)/n))
      Je[i]=sum(r[-i]*wi) 
    }
    Jve=(1-mean(pi))*((n-1)/n)*sum((Je-mean(Je))^2)
  }
  
  if(((strata==FALSE)&(cluster==TRUE))|((strata==TRUE)&(cluster==TRUE))){
    Je=vector()
    for(i in levels(clu)){
     wi=1/(pi[clu!=i]*((n-1)/n))
     Je[i]=sum(r[clu!=i]*wi) 
    }
    Jve=(1-mean(pi))*((n-1)/n)*sum((Je-mean(Je))^2)
  }

  if(Jve<0){warning("The variance estimation can not be negative.")}
  return(Jve)
}
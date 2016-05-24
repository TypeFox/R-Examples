# Package Code 'fishmove'
# 
# Author: Johannes Radinger
###############################################################################


#Function for queryinf fish dispersal kernel (reverse of fishmove)
fishmove.query <-function(fishmove,p=0.67,dist=NA,fromto=NA,prob=NA,reach=NA,w=1,level="fit",...){
  
    #check of dimensions of fishmove
    if(length(fishmove$pred.fishmove["fit","sigma_stat",,,,])>1 && !hasArg(dim)){
      warning("Multiple values supplied in fishmove. Only first values used.",call. = FALSE)
    }
    
    fishmove.array <- fishmove$pred.fishmove[,,1,1,1,1]
    
    f <- function(dist,
                  sigma_stat=fishmove.array[level,"sigma_stat"],
                  sigma_mob=fishmove.array[level,"sigma_mob"]) {
      (pnorm(q=dist,mean=0,sd=sigma_stat)*p+pnorm(q=dist,mean=0,sd=sigma_mob)*(1-p))*w
    }
    
    # if probability at a given dist is queried
    if(!is.na(dist)){
      out <- f(dist=(abs(dist))*-1)
    }
    
    #if probability between two distancess is queried 
    if(!is.na(fromto)){
      if(!is.vector(fromto) && length(fromto)==2){stop("Wrong format for fromto provided. Needs vector of length = 2")}
      out <- abs(f(dist=fromto[1])- f(dist=fromto[2])) # is abs here the correct way?
    }
    
    #if distance is queried based on a probability
    if(!is.na(prob) && is.na(reach)){
    out <- uniroot(function(x) f(dist=x)-prob,c(-1e25,0))$root
    }
    
    #if distance is queried based on a probability and reach
    if(!is.na(prob) && !is.na(reach)){
      if(prob>f(dist=reach)-f(dist=0)) {stop(paste("The provided prob",prob,"is higher than possible for a reach of",reach,"m length.",sep=" "))}
      else out <- uniroot(function(x) f(dist=x+(reach/2))-f(dist=x-(reach/2))-prob,c(-1e25,0))$root
    }
       
    return(out)
    
}


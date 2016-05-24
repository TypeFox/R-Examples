getSubclasses <- function(className){
  class <- getClass(className)
    
  subclasses <- names(class@subclasses)
  index <- sapply(subclasses,function(x){!getClass(x)@virtual})
  
  names <- subclasses[index]  
  
  if(!class@virtual){
    names <- c(className,names)    
  }
  
  return(names)
}

mise <- function(model1,model2,discreteApproximation = TRUE){
  
  if(is(model1,"BoundedDensity") & is(model2,"BoundedDensity")){
  
    if(discreteApproximation){
    
      if(!any(getdataPointsCache(model1) != getdataPointsCache(model2))){
        sum((getdensityCache(model1)-getdensityCache(model2))^2)/length(getdensityCache(model1))
      }
      else{
        stop("DataPointsCache must be the same for both models when discreteApproximation = TRUE. Use discreateApproximation == FALSE instead")
      }
      
    }
    else{ # no discrete approximation
      dif <- function(x,model1,model2){
        (density(model1,x) - density(model2,x))^2
      }
      integrate(dif,lower=0,upper=1,model1=model1,model2=model2)
    }
        
  }
  else{
    stop("model1 and model2 must be BoundedDensity objects. See getSumclasses(\"BoundedDensity\") to see all the allowed classes")
  }
  
}
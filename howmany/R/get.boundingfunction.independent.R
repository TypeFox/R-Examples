"get.boundingfunction.independent" <-
function(m,alpha,at=(1:1000)/1000,method="asymptotic")
  {
    if(method!="asymptotic") stop(paste("Method ", method, " not available"))
    beta <- find.beta(m,alpha)

    if( max(at)>1 | min(at)<0 ) stop( "bounding-function can only be evaluated within [0,1]")
    if(length(at)<1) stop(" bounding-function must be evaluated at least at one point ")
    
    boundingfunction <-   m*(  at + beta*sqrt(at*(1-at))  ) 
    return(boundingfunction)
        
  }


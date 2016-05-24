`checkCenter` <- function(center,X,...){
   k <- ncol(X) 
   if(is.numeric(center))
   {
     if(length(center)!= k) stop("'center' has a different dimension than 'X'")
     CENTER <- center
   }
 
   if(is.function(center))
   {
     CENTER <- do.call(center, list(X, ...))
     if(!is.vector(CENTER) | length(CENTER)!=k) 
       stop("'center' does not return a vector with the same dimension as 'X'")
   }

   if(is.character(center))
   {
     center <- match.arg(center,c("colMean", "ojaMedian", "spatialMedian", "compMedian", "HRMedian"))
     CENTER<-switch(center,
        "colMean"={colMeans(X,...)}
        ,
        "ojaMedian"={ojaMedian(X,...)}
        ,
        "spatialMedian"={spatial.median(X, ...)}
        ,
        "compMedian" = {apply(X,2,median,...)}
        ,
        "HRMedian"={HR.Mest(X, ...)$center}
        )
   }
   return(CENTER)
}

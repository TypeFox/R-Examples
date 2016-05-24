#### Distance Metrics

#' @export

distance <- function(x,y,method="euclidean",tvectors=tvectors,breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){
    
    if(breakdown==TRUE){
      
      x <- breakdown(x)
      
      y <- breakdown(y)
      
    }
    
    
    if(x %in% rownames(tvectors) && y %in% rownames(tvectors)){
      
      v <- tvectors[x,]
      w <- tvectors[y,]
      
      if(method == "euclidean"){
        
        dist <- sqrt( sum( (v - w)^2 ) )
        
      }
      
      
      if(method == "cityblock"){
        
        dist <- sum( abs(v - w))
        
      }
      
      
    }else{
      
      warning("x and y must be words in rownames(tvectors)")
      return(NA)
      
    }
    
  }else{stop("tvectors must be a matrix!")}
  
}
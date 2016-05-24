#### Easier Two-Word Cosine 

#' @export
#' @importFrom lsa cosine 

Cosine <- function(x,y,tvectors=tvectors,breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){
    
    if(breakdown==TRUE){
      
      x <- breakdown(x)
      
      y <- breakdown(y)
      
    }
    
    
    if(x %in% rownames(tvectors) && y %in% rownames(tvectors)){
      
      as.numeric(cosine(tvectors[x,],tvectors[y,]))
      
    }else{
      
      warning("x and y must be words in rownames(tvectors)")
      return(NA)
      
    }
    
  }else{stop("tvectors must be a matrix!")}
  
}
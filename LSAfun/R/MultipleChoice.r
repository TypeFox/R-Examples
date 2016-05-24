
#' @export
MultipleChoice <- function(x,y,tvectors=tvectors,breakdown=FALSE){

 
  cos <- vector(length=length(y))
  names(cos) <- y
  
  for(j in 1:length(y)){value <- costring(x,y[j],
                                           tvectors=tvectors,
                                           breakdown=breakdown)
  
  if(is.na(value)){cos[j] <- -2}
  if(!is.na(value)){cos[j] <- value}                      
                        
  }
  
  
  names(cos)[which(cos == max(cos))]
}

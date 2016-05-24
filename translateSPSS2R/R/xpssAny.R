#' Selecting cases by condition
#'
#' xpssAny can be perceived as a wrapper function for \code{\link{\%in\%}} applicable on more than one variable.
#'
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @param st atomic numeric or atomic character with a single value to search for OR variable where to search in. 
#' @param nd atomic numeric or atomic character, respecetively as numeric vector or character vector with values to search for OR variables to search in. 
#' @return A logical vector with matched conditions.
#' @author Andreas Wygrabek
#' @seealso \code{\link{xpssCount}} \code{\link{\%in\%}} \code{\link{is.element}}
#' @examples
#' 
#' data(fromXPSS)
#' 
#' xpssAny(fromXPSS, 310, c("V7_1", "V7_2"))
#' 
#' xpssAny(fromXPSS, "V7_1", c(310,320,170))
#'
#' xpssAny(fromXPSS, "Audi", c("V1", "V7_2"))
#' 
#' xpssAny(fromXPSS, "V1", c("Audi"))
#' @export
xpssAny <- function(x, st = NULL, nd = NULL) {
    
    functiontype  <- "SB"
    x <- applyMetaCheck(x)
  
    temp <- data.frame(1:length(x=x[,1]))
    
    if(length(st)>1)
    {
      for(j in 1:length(st)) {
        if(!(st[[j]] %in% colnames(x))){
          
          logDF <- matrix(nrow = nrow(x), ncol = length(nd))    
          
          for(i in 1:length(nd)){
            logDF[,i] <- x[,nd[i]] %in% st[[j]] 
          }
          
          temp[[j]] <- rowSums(logDF) > 0
          
        } else { 
          
          temp[[j]] <- x[,st[[j]]] %in% nd
          
        }
      }
         out <- rowSums(temp)
        for(i in 1:length(out)) {
          if(out[[i]]>0){
            out[[i]] <- TRUE
          } else {
            out[[i]] <- FALSE
          }
        }       
       out <- as.logical(out)
    } else {
      if(!(st %in% colnames(x))){
        
        logDF <- matrix(nrow = nrow(x), ncol = length(nd))    
        
        for(i in 1:length(nd)){
          logDF[,i] <- x[,nd[i]] %in% st 
        }
        
        out <- rowSums(logDF) > 0
        
      } else { 
        
        out <- x[,st] %in% nd
        
      }
    }
    return(out)
}





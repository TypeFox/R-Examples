#' Merges Attributdata 
#'
#' Attribut merge function
#'
#' Helper Function for DeMerging Attributedatasets
#' 
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @return Output is the input data demerged on the basis of the attributes of the actual dataset.
#' @author Bastian Wiessner
#' @keywords internal

applyAttributeDemerge <- function(x) {
  
    if(!is.null(attributes(x)$FILTER)){ 
      if(attributes(x)$FILTER != FALSE && "DM" %in% class(x) == T) {
          
        pos <-  which(eval(parse(text=paste("x$",attributes(x)$FILTER,sep=""))))
        
                
        attr(x, "FILTERED_DATA") <- x[setdiff(1:nrow(x),pos),]
        
        attBack <- attributesBackup(x)
        
        x <- x[pos,]

        x <- applyAttributes(x, attBack)
      }
    } 
    
  return(x)
}


#' Apply stored attributes 
#'
#' Applies attributes stored by attributesBackup
#'
#' @usage applyAttributes(x, attributesToApply = NULL)
#' @param x a data.frame or input data of class xpssFrame. 
#' @param attributesToApply to applied attributes
#' @return Object with attributes from \code{\link{attributesBackup}}
#' @seealso \code{\link{attributes}} \code{\link{attr}}
#' @author Andreas Wygrabek
#' @examples
#' #load data
#' data(fromXPSS)
#' attributes(fromXPSS)
#' attributes(fromXPSS$V7_2)
#' x <- attributesBackup(fromXPSS)
#' fromXPSS <- fromXPSS[order(fromXPSS$V2),] 
#' attributes(fromXPSS)
#' attributes(fromXPSS$V7_2)
#' fromXPSS <- applyAttributes(fromXPSS, x)
#' attributes(fromXPSS)
#' attributes(fromXPSS$V7_2)
#' @export


applyAttributes <- function(x, attributesToApply = NULL){
    
    
    if(!"xpssAttributes" %in% class(attributesToApply))
    {
        stop("Attributes have to be stored in an object of class xpssAttributes")
    }
    

    
    ####
    
    backRN <- rownames(x)

    attributesToApply[["global"]]$row.names <- backRN
    #----
    
    attributes(x) <- attributesToApply[["global"]] 
    

    for(i in 1:ncol(x)){
        if((!is.null(attributes(x[,i])$varname)) && !is.null(names(attributesToApply$local))){
        attributes(x[,i]) <- attributesToApply$local[[which(names(attributesToApply$local) == attributes(x[,i])$varname)]]
        } else if(names(x)[i] %in% names(attributesToApply$local)){
            attributes(x[,i]) <- attributesToApply$local[[which(names(attributesToApply$local) %in% names(x)[i])]] 
            } else if(colnames(x)[i] %in% colnames(attributesToApply$local)){
              attributes(x[,i]) <- attributesToApply$local[,which(colnames(attributesToApply$local) %in% colnames(x)[i])]
            }
    }
    
    return(x)
    
}

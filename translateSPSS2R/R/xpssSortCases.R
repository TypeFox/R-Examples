#' Sorts data ascending or descending order
#'
#' R implementation of the SPSS \code{SORT CASES} argument. xpssSortCases reorders the sequence of cases in the dataset based on the values of one or more variables. 
#' 
#' The argument order has to be of the same length as the argument variables. Optionally, the sorting can be specified in ascending or descending order for any variable. It is also possible to use combinations of ascending and descending order for different variables.
#'
#' @usage xpssSortCases(x, variables = NULL, order ="A")
#' @param x a (non-empty) data.frame or input data of class "xpssFrame".  
#' @param variables atomic character or character vector with the names of the variables. Also \code{\link{rownames}} can be used to sort the data. 
#' @param order atomic character or character vector containing either "A" for ascending order or "D" for descending order. 
#' @return Returns a sorted xpssFrame.
#' @author Andreas Wygrabek
#' @seealso \code{\link{sort}} \code{\link{order}}
#' @examples
#' 
#' data(fromXPSS)
#' 
#' xpssSortCases(fromXPSS, variables = c("V4", "V7_1", "V7_2"), order = c("A","D","A"))
#' @export
xpssSortCases <- function(x, variables = NULL, order = "A"){
   
    stopifnot(length(variables) == length(order))  
    stopifnot(order == "A" | order == "D" | order == "UP" | order == "DOWN")
    
    ####################################################################
    ####################################################################
    
    functiontype <- "DM"
    x <- applyMetaCheck(x)
    sortframe <- data.frame(x)
    ####################################################################
    ####################################################################
    ####################################################################
    
    attBack <- attributesBackup(x)
    
    if("rownames" %in% variables){
        
        x <- x[order(rownames(x)),]
    } else {
      
      #Do: Sort Ascending If order is 'A' or 'UP', else sort descending
      for(i in 1:length(variables)){
        if(order[[i]] == "D" && class(x[,variables[i]]) == "character"){
          eval(parse(text=paste(paste0("sortframe$",variables[i])," <- -as.numeric(as.factor(x[,",which(names(x) %in% variables[i]),"]))")))
        } else{
          eval(parse(text=paste(paste0("sortframe$",variables[i])," <- x[,",which(names(x) %in% variables[i]),"]")))
        }        
      }
      
    eval(parse(text = paste("x <- x[order(",paste0("sortframe[,",1:length(order),"]",collapse=","),"),]")))
    
    }
    x <- applyAttributes(x, attBack)
    
    x <- applyAttributeDemerge(x)
        
    return(x)   
}

#-----------------------------------------------------------------------------#
#                                                                             #
#            QUALITY CONTROL AND RELIABILITY IN R                             #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sánchez                                       #
#              Student Master of Statistical Techniques                       #
#              University of The Coruña, SPAIN                                #
#              mflores@outlook.com                                            #
#                                                                             #
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Main function to create a 'mqcd' object
#-----------------------------------------------------------------------------#
##' Multivariante Quality Control Data
##' 
##' Create an object of class 'mqcd' to perform statistical quality control.
##' This object be used to plot Multivariate Control Charts.
##' 
##' @aliases mqcd 
##' @param data a matrix or data-frame or array where it should contain data.
##' @param data.name  a string that specifies the title displayed on the plots. 
##' If not provided is taken from the name of the object data.
##' @export
##' @examples
##' library(qcr)
##' data(dowel1)
##' str(dowel1)
##' data.mqcd <- mqcd(dowel1)
##' str(data.mqcd)


mqcd <- function(data, data.name = NULL)
  #.........................................................................
{

  if (!is.matrix(data) & !is.data.frame(data) & !is.array(data))
    stop("object must be a matrix or data.frame or array")
  
  p <- ncol(data) # quality characteristics
  m <- nrow(data) # number of samples or observations
  if (class(data) == "matrix" || class(data) == "data.frame") {
      names <- colnames(data)    
      data <- array(data.matrix(data),c(m,p,1))
      colnames(data) <- names        
  }    
  n <- dim(data)[3] # observations or sample size 
  
  if (is.null(data.name))
    data.name <- " DATA"
  
  result <- data
  
  attr(result, "data.name") <- data.name
  attr(result, "type.data") <- "Multivariate"
  
  oldClass(result) <- c("mqcd", "array")
  
  
  return(result)
} # mqcd
#.........................................................................

################################################################################
#  Print method for class 'select.parfm'                                       #
################################################################################
#                                                                              #
#  This function prints the objects of class 'parfm'                           #
#                                                                              #
#  Its parameters are                                                          #
#   - x         : the object of class 'select.parfm'                           #
#   - digits    : number of significant digits                                 #
#   - na.prints : character string indicating NA values in printed output      #
#                                                                              #
#                                                                              #
#   Date:                 January, 10, 2012                                    #
#   Last modification on: June 27, 2012                                        #
################################################################################

print.select.parfm <- function(x,
                               digits=3,
                               na.print="----",
                               ...) {
  if (missing(digits))
    digits <- 3
  
  if (!is.null(x)){
    cat("\nAIC:\n")  
    print(round(x$AIC, digits), na.print=na.print)
    
    cat("\n\nBIC:\n")  
    print(round(x$BIC, digits), na.print=na.print)
  } 
}

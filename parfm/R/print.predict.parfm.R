################################################################################
#  Print method for class 'predict.parfm'                                      #
################################################################################
#                                                                              #
#  This function prints the objects of class 'predict.parfm'                   #
#                                                                              #
#  Its parameters are                                                          #
#   - x         : the object of class 'select.parfm'                           #
#   - digits    : number of significant digits                                 #
#   - na.prints : character string indicating NA values in printed output      #
#                                                                              #
#                                                                              #
#   Date: January, 10, 2012                                                    #
#                                                                              #
################################################################################

print.predict.parfm <- function(x,
                                digits=3,
                                na.print="",
                                ...) {
  if (!is.null(x)){
    frailty <- attr(x, "frailty")
    dist <- attr(x, "dist")
    cat(paste(toupper(substr(frailty, 1, 1)), substr(frailty, 2, 100), 
                    " frailty model with ", 
                    toupper(substr(dist, 1, 1)), substr(dist, 2, 100),
                    " baseline\n", sep="", collapse=""))
    
    toprint <- as.matrix(cbind(names(x), round(as.vector(x), digits)))
      colnames(toprint) <- c(attr(x, "clustname"), "frailty")
      rownames(toprint) <- rep("", nrow(toprint))
    print(toprint, quote=FALSE)
  } 
}

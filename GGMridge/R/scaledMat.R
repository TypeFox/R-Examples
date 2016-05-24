#----------------------------------------------------------------------#
# Scale a square matrix to have unit diagonal elements.                #
#----------------------------------------------------------------------#
#                                                                      #
# Inputs :                                                             #
#                                                                      #
#  x A square matrix with positive diagonal elements                   #
#                                                                      #
# Outputs :                                                            #
#                                                                      #
# The scaled matrix                                                    #
#                                                                      #
#----------------------------------------------------------------------#
scaledMat <- function(x) {

  diagx <- diag(x)

  newx <- x / sqrt(diagx %o% diagx)

  return(newx)
}

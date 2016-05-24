
###############################################################################
# CheckInput:                                                                #
#                                                                             #
# Check the consistency of the input arguments and whether the format meets   #
# requirements.                                                               #
###############################################################################

CheckInput <- function(response, Q) {
  
  nperson <- nrow(response)
  nitem <- ncol(response)
  nitem.Q <- nrow(Q)
  
  if (nitem != nitem.Q) {
    return("Item numbers in the response matrix and Q-matrix do not agree.")
  }
  
  if (!all(response %in% c(1, 0))) {
    return("The response matrix should have only two values: 1=correct, 0=incorrect.")
  }
  
  if (!all(Q %in% c(1, 0))) {
    return("The Q-matrix should have only two values: 1=attribute is required, 0=attribute is not required.")
  }
    
}

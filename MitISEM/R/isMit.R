# Function to check if mixture of multivariate t distributions is well-defined
#
# inputs: 
#   mit : [list] of distribution parameters
# outputs:
#   logical indicating success (TRUE) or failure (FALSE)
#
# author : Nalan Basturk
# date   : 20120912

isMit<-function(mit){
  # list/vector/matrix input/components
  r1    <- is.list(mit) 
  if(r1)
    r1    <- (is.vector(mit$p)&is.matrix(mit$mu)&is.matrix(mit$Sigma)&is.vector(mit$df))
  # vector/matrix dimensions
  if(r1){
    H.mit <- length(mit$p)
    k.mit <- ncol(mit$mu)
    r1    <- (H.mit>0 & length(mit$p)==H.mit & length(mit$df)==H.mit & nrow(mit$mu)==H.mit & nrow(mit$Sigma)==H.mit)
    r2    <- (ncol(mit$mu)==k.mit & ncol(mit$Sigma==k.mit))    
    r1    =  (r1 & r2)
  }
  as.logical(r1)
}

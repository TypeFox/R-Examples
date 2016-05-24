#-----------------------------------------------------------------------------------------------------------------#
# Function modelLSVM
#-----------------------------------------------------------------------------------------------------------------#

#Return the parameters of the hyperplanes used in the construction of the LSVM classifier
#' @export
modelLSVM <- function(X, Y, convexity){

  if( convexity == 0){
    stop('The problem is not convex !')
  }

  if( abs(convexity) != 1){
    stop('convexity must be equal to -1 or +1 !')
  }  


  XS <- X[which(Y > 0), ]
  XF <- X[which(Y <= 0), ]
  d <- dim(X)[2]

  if(is.null(dim(XF))){
    rXF <- 1
    XF <- matrix(XF, nrow = 1)
  }else{
    rXF <- dim(XF)[1]
  }

  if(is.null(dim(XS))){
    rXS <- 1
    XS <- matrix(XS, nrow = 1)
  }else{
    rXS <- dim(XS)[1]
  }

  #Build a set of hyperplane separating each point of XS from all the point of XF
  if(convexity == -1){
    A <- matrix(0, ncol = d + 1, nrow = rXS)
    for(i in 1:rXS){
      A[i, ] <- SVMLinearMonotone(rbind(XS[i,], XF), c(1, rep(-1, rXF))  )
    }
  }

  #Build a set of hyperplane separating each point of XF from all the point of XS
  if(convexity == +1){
    A   <- matrix(0, ncol = d + 1, nrow = rXF)
    for(i in 1:rXF){
      A[i, ] <- SVMLinearMonotone(rbind(XF[i,], XS), c(-1, rep(1, rXS))  )
    }
  }

  return(A)
}


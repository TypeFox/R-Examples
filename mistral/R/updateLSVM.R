#-----------------------------------------------------------------------------------------------------------------#
# Function updateLSVM
#-----------------------------------------------------------------------------------------------------------------#

#Update the classifier
#' @import emoa
#' @export
updateLSVM <- function(X.new, Y.new, X, Y, A.model.lsvm, convexity, PLOTSVM = FALSE, step.plot.LSVM = 1, hyperplanes = FALSE, limit.state.estimate = TRUE){

  if( sum(sign(Y.new) == convexity) > 0 ){
    return( modelLSVM( rbind(X, X.new), c(Y, Y.new), convexity))
  }
  
  d <- dim(X)[2]

  if(is.null(dim(X.new)) ){
    rX.new <- 1
    X.new <- matrix(X.new, nrow = 1)
  }else{
    rX.new <- dim(X.new)[1]
  }
  #Add hyperplanes, built with the new data, to A.model.lsvm 
  if(convexity == -1){
    A   <- matrix(0, ncol = d + 1, nrow = rX.new)
    XF  <- X[which(Y <= 0), ]
    rXF <- ifelse( is.null(dim(XF)), 1, dim(XF)[1] )
    for(i in 1:rX.new){
      A[i, ] <- SVMLinearMonotone(rbind(X.new[i,], XF), c(1, rep(-1, rXF))  )

      if(PLOTSVM){
        if( (i %% step.plot.LSVM) == 0){
          plotLSVM(X, Y, A.model.lsvm, hyperplanes = FALSE, limit.state.estimate = TRUE, convexity = convexity)
          points(X.new[i, 1], X.new[i, 2], pch=19, col='green')
          continue <- readline("Continue?(y/n) ")
          if(continue != 'y'){
            PLOTSVM <- FALSE
          }
        }
      }
    }
  }

  if(convexity == +1){
    A   <- matrix(0, ncol = d + 1, nrow = rX.new)
    XS  <- X[which(Y > 0), ]
    rXS <- ifelse( is.null(dim(XS)), 1, dim(XS)[1] )
    for(i in 1:rX.new){
      A[i, ] <- SVMLinearMonotone(rbind(X.new[i,], XS), c(-1, rep(1, rXS))  )
      if(PLOTSVM){
        if( (i %% step.plot.LSVM) == 0){
          plotLSVM(X, Y, A.model.lsvm, hyperplanes = FALSE, limit.state.estimate = TRUE, convexity = convexity)
          points(X.new[i,1], X.new[i,2], pch=19, col='green')
          continue <- readline("Continue?(y/n)")
          if(continue != 'y'){
            PLOTSVM <- FALSE
          }
        }
      }
    }
  }
  return(rbind(A.model.lsvm, A))
}


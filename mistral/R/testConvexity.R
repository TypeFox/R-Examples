#-----------------------------------------------------------------------------------------------------------------#
# Function testConvexity
#-----------------------------------------------------------------------------------------------------------------#

# Test if one of the class is convex and proiveds the hyperplanes used to build the LSVM
#' @export
testConvexity <- function(X,Y){

  XS <- X[which(Y > 0), ]
  XF <- X[which(Y <= 0), ]
  d  <- dim(X)[2]

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


  A.low <- matrix(0, ncol = d + 1, nrow = rXS)
  test.low <- TRUE
  i <- 1
  #test if each points of XS can be separated from all points of XF by a hyperplane
  while(test.low){
    res.low <- SVMLinearMonotone(rbind(XS[i,], XF), c(1, rep(-1, rXF)))
    A.low[i, ] <- res.low
    if( length(res.low) == 1){
      test.low <- FALSE
    }
    i <- i + 1
    if(i > rXS){
      cat('The set associated to the label -1 is convexe \n')
      return( list(-1, A.low) )
    }
  }

  #test if each points of XF can be separated from all points of XS by a hyperplane
  A.up <- matrix(0, ncol = d + 1, nrow = rXF)
  test.up <- TRUE
  i <- 1
  while(test.up){
    res.up    <- SVMLinearMonotone(rbind(XF[i,], XS), c(-1, rep(1, rXS))  )
    A.up[i, ] <- res.up
    if( length(res.up) == 1){
      test.up <- FALSE
    }
    i <- i + 1
    if(i > rXF){
      test.up <- FALSE
      cat('The set associated to the label -1 is concave \n')
      return(list(1, A.up))
    }
  }

  cat('The problem is not convex!')
  return(0)
}




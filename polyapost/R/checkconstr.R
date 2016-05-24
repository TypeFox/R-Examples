#checks if the matrices and vectors defining the
#constaints are set up properly

checkconstr<-function(A1,A2,A3,b1,b2,b3)
  {
    if( !is.matrix(A1)) stop("A1 is not a matrix")
    if( !is.vector(b1)) stop("b1 is not a vector")
    if( min(b1)< 0) stop("b1 has negative values")
    if( nrow(A1) != length(b1)) stop("nrow(A1) != length(b1)")
    if( !is.matrix(A2)) stop("A2 is not a matrix")
    if( !is.vector(b2)) stop("b2 is not a vector")
    if( min(b2)< 0) stop( "b2 has negative values")
    if( nrow(A2) != length(b2)) stop("nrow(A2) != length(b2)")
    if( !is.null(A3)){
      if( !is.matrix(A3)) stop("A3 is not a matrix")
      if( !is.vector(b3)) stop("b3 is not a vector")
      if( min(b3)< 0) stop("b3 has negative values")
      if( nrow(A3) != length(b3)) stop("nrow(A3) != length(b3)")
    }
  }

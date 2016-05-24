matexp <- function(a, dt=1, order=8) {
    if (!is.matrix(a))
        stop("invalid (non-matrix) argument")
    da <- dim(a)
    if (da[1] != da[2])
      stop("matrix not square")
    if (!is.numeric(order) | order < 1 )
      stop("order must be a positive integer number")
    
    # What if the matrix has only zero elements or dt=0
    if(dt==0 | all(a == 0) ) 
      return(diag(da[1]))
    
    # Internals
    # SUBROUTINE DGPADM( IDEG,M,T,H,IFLAG )
    Fobj <- .Fortran("MATEXPFORTRANSUB",
                     as.integer(order), #IDEG  1
                     as.integer(da[1]), #M     2
                     as.double(dt),     #T     3
                     as.double(a),      #H     4
                     as.integer(0),     #IFLAG 5
                     PACKAGE="PSM")
    if(Fobj[[5]] < 0)
      stop("Unable to determine matrix exponential")
    return( matrix(Fobj[[4]], nrow=da[1]))
}



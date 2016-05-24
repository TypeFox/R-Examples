L2Distance <- function(a,b,df = 0){
#    A - (DxM) matrix 
#    B - (DxN) matrix
#    df = 1, force diagonals to be zero; 0 (default), do not force
# 
# Returns:
#    E - (MxN) Euclidean distances between vectors in A and B
#
#
# Description : 
#    This fully vectorized (VERY FAST!) m-file computes the 
#    Euclidean distance between two vectors by:
#
#                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
#
# Example : 
#    A = rand(400,100); B = rand(400,200);
#    d = distance(A,B);

# Author   : Roland Bunschoten
#            University of Amsterdam
#            Intelligent Autonomous Systems (IAS) group
#            Kruislaan 403  1098 SJ Amsterdam
#            tel.(+31)20-5257524
#            bunschot@wins.uva.nl
# Last Rev : Wed Oct 20 08:58:08 MET DST 1999
# Tested   : PC Matlab v5.2 and Solaris Matlab v5.3

# Copyright notice: You are free to modify, extend and distribute 
#    this code granted that the author of the original code is 
#    mentioned as the original author of the code.

# Fixed by JBT (3/18/00) to work for 1-dimensional vectors
# and to warn for imaginary numbers.  Also ensures that 
# output is all real, and allows the option of forcing diagonals to
# be zero.  

  ## if a is a vector convert it to a matrix
  if(is.matrix(a)==FALSE){
    a <- t(matrix(a))
  }
  ## if b is a vector convert it to a matrix
  if(is.matrix(b)==FALSE){
    b <- t(matrix(b))
  }


  if(dim(a)[1] != dim(b)[1]){
    stop("L2distance: A and B should be of same dimensionality")
  }

  
  if (all(is.numeric(a), is.numeric(b)) == FALSE){
    stop('Warning: running distance.m with imaginary numbers.  Results may be off.')
  }

  D <- dim(a)[1]
  if(D == 1){
    a <- rbind(a,rep(0,D))
    b <- rbind(b,rep(0,D))
  }

  aa <- matrix(apply(a*a,2,sum)); bb <- matrix(apply(b*b,2,sum)); ab <- t(a)%*%b


  ## copied from: http://haky-functions.blogspot.com/search/label/matlab2R
  repmat = function(X,m,n){
    ##R equivalent of repmat (matlab)
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=TRUE)
  }


  ## look that elements are either zero or positive
  d <- sqrt(apply(repmat(aa,1,dim(bb)[1]) + repmat(t(bb),dim(aa)[1],1)-2*ab,c(1,2),function(x){if(x>0){x}else{0}}))

  
  # force 0 on the diagonal? 
  if (df==1){
    diag(d) <- 0
  }
  return(d)
}

`AXB.prod` <-
function( A,  B)
  {
#######   cross product of A cross B in 3D
    ##  A and B can be lists of x,y,z or vectors
    ####   returns a vector

    ###  example:
  ####   AXB.prod(c(1,0,0), c(0,1,0)) =  c(0, 0 , 1)
 ####    AXB.prod(c(0,1,0), c(1,0,0)) =  c(0, 0 ,-1)


    if(is.list(A)) A = unlist(A)
    if(is.list(B)) B = unlist(B)

    if(length(A)<3) A = c(A, 0)
    if(length(B)<3) B = c(B, 0)
    
    c = rep(0, 3)
    
    c[1] = (A[2]*B[3]) - (B[2]*A[3]);
    c[2] = (B[1]*A[3]) - (A[1]*B[3]);
    c[3] = (A[1]*B[2]) - (B[1]*A[2]);
    return (c);
  }


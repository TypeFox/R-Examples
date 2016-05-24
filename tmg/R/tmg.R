rtmg <- function(n, M, r, initial, f = NULL, g=NULL, q=NULL, burn.in=30){

  d = nrow(M)      # dimension of target space
  if (ncol(M)!=d){
    cat("Error: M must be a square matrix.")
    return()
  }  
  
  if (length(initial)!=d){
    cat("Error: wrong length for initial value vector.")
    return()
  }
    
  # symmetrize M and verify that it is positive definite  
  M = (M + t(M))/2;  
  eigs=eigen(M, symmetric=T, only.values=T)$values
  if (any(eigs<=0)){
    cat("Error: M must be positive definite.")
    return()
  }
    
  
  # we change variables to the canonical frame, sample by calling the C++ code and transform back to the original frame

  R = chol(M)
  Mir = solve(M,r)
  Ri = solve(R)
  initial2= as.vector(R%*%initial) - as.vector(r%*%Ri) 
  
  if (is.null(f) & is.null(g)){
    numlin=0
    f2=NULL
    g2=NULL    
    } else   # if there are linear constraints
      { if (is.matrix(f) & is.vector(g)){  
    # verify linear constraints sizes    
    numlin = nrow(f);
    if (length(g) != numlin | ncol(f) != d ){
      cat("Error: inconsistent linear constraints. f must be an m-by-d matrix and g an m-dimensional vector.")
      return()
    }
  
    # verify initial value satisfies linear constraints
    if (any(f%*%initial + g <= 0)) {
      cat("Error: initial point violates linear constraints.")
      return()
    }
    
    # map linear constraints to the canonical frame
    f2 = f%*%Ri
    g2 = as.vector(f%*%Mir + g)
    
  } # if (is.matrix(f) & is.vector(g))
     else { 
    cat("Error: for linear constraints, f must be a matrix and g a vector.\n")
    return()
  }} 
  
  
  if (is.null(q)){
    numquad=0 
    quads = NULL
    }
  else    # if there are quadratic constraints
      { if (is.list(q)){
        # verify that the elements of the quadratic constraints are lists of length  3
        ll=lapply(q,length)
        nc =c(do.call("cbind",ll))
        if (any(nc!=3)) {
          cat("Error: each element in q must be a list of length 3.\n");
          return()
        }
  
        numquad = length(q);
        quads = matrix( nrow=numquad*(d+2), ncol=d)
     
        for (i in 1:numquad){
          qci = q[[i]];
          t = initial %*% qci[[1]] %*% initial + qci[[2]] %*% initial + qci[[3]];
     
          if (t <=0){
            cat("Error: initial point violates quadratic constraints. \n")
            return()
          } else {
            
            # map quadratic constraints to the canonical frame
            A= qci[[1]]
            B= qci[[2]]
            C= qci[[3]]
            quads[ ((i-1)*(d+2)+1): ((i-1)*(d+2)+d)  , ] = t(Ri)%*% A %*% Ri
            quads[ ((i-1)*(d+2)+d+1) ,   ] = 2*t(Mir) %*% A %*% Ri + t(B) %*% Ri 
            C = C + t(Mir)%*%A%*% Mir + t(B)%*% Mir
            quads[ i*(d+2) ,  ] = c(C, rep(0,d-1))          
          }        
      } #for (i in 1:numquad)       
     } #if (is.list(q))
     else {
       cat("Error: for quadratic constraints, q must be a list.\n")
       return()
  }}
  
  seed = sample(1:10000000,1)

  samples = .Call( "rtmg", n+burn.in, seed,  initial2, numlin, f2, g2, numquad, quads, PACKAGE = "tmg" )
  samples = samples[(burn.in+1):(burn.in+n),]
  samples = samples%*%t(Ri) + matrix(rep(Mir,n),nrow=n, byrow=T)
  return(samples)
}

getw <-
function(I,Y, Isign= -1*(abs(I)>10^(-5)) ,wleafs=NULL, rho = 10^(-5), epsilon=Inf, silent=FALSE ){

  Dmat <- t(I) %*% I

  pI <- ncol(I)
  nI <- nrow(I)
  
    
  sv <- eigen( t(Isign)%*%Isign )
  rho <- max(c( rho, -min(sv$val) + rho))
  ind <- which( sv$val < 10^(-5) )
  c <- - as.numeric( t( Y - I%*%wleafs ) %*% I ) 
  
  if( length(ind) > 0){
    if(!silent) cat("\n dimension of null space of I                           :",length(ind))
    basis <- sv$vec[,ind,drop=FALSE]
    csmall <- as.numeric(c%*%basis) 
    Hsmall <- t(basis) %*% Dmat %*% basis
    Hsmall <- 0.5*( Hsmall +t(Hsmall))

    if(epsilon==Inf){
      sqp <- solve.QP(Hsmall + rho * diag(nrow(Hsmall))  , -csmall, t(basis), - 0.99*wleafs - 10^(-8)  )
    }else{
      if(notpenroot <- FALSE) epsconstraint <- apply( basis[1:(pI-1),],2,sum)
      epsconstraint <- apply( basis,2,sum)
      sqp <- solve.QP(Hsmall + rho * diag(nrow(Hsmall))  , -csmall, t(rbind( basis, -epsconstraint)) , c( -wleafs - 10^(-7), -epsilon)  )
           
    }
    
    w <-  wleafs +  basis %*% sqp$solution

    
  }else{
    warning("trivial null space -- increase number of nodes to get meaningful results")
    w <- wleafs
  }
  if(!silent) cat("\n number of selected nodes                               :",sum(w>0.01),"\n")

  
  return(w)
  
}


getw2 <-
function(I,Y, Isign= -1*(abs(I)>10^(-5)) ,wleafs=NULL, rho = 10^(-5), epsilon=Inf, silent=FALSE ){

  Dmat <- t(I) %*% I
  Ismat <- t(Isign)%*%Isign

  pI <- ncol(I)
  nI <- nrow(I)
  
   
  c <- - as.numeric( t( Y - I%*%wleafs ) %*% I ) 
   

  delta <- solve.QP( Dmat + max(abs(I)) *20 * Ismat + diag(pI), -c,  t( rbind( rep(1,pI), rep(-1,pI), diag(pI))), c(-0.001,-0.001,-wleafs))$solution
    
    
  w <-  wleafs +  delta

    
  if(!silent) cat("\n number of selected nodes                               :",sum(w>0.01),"\n")

  
  return(w)
  
}


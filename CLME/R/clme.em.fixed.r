#' @description \code{clme_em_fixed} performs a constrained EM algorithm for linear fixed effects models.
#' 
#' @rdname clme_em
#' @export
#' 
clme_em_fixed <- function( Y, X1, X2 = NULL, U = NULL, Nks = dim(X1)[1],
                     Qs = dim(U)[2], constraints, mq.phi = NULL, tsf = lrt.stat, 
                     tsf.ind = w.stat.ind, mySolver="LS", em.iter = 500, 
                     em.eps =  0.0001, verbose = FALSE, ... ){
  
  if( verbose==TRUE ){
    message("Starting EM algorithm")    
  }
  
  N  <- sum(Nks)
  N1 <- 1 + cumsum(Nks) - Nks
  N2 <- cumsum(Nks)
    
  X  <- as.matrix(cbind(X1, X2))
  theta.names <- NULL
  if( is.null(colnames(X))==FALSE ){
    theta.names <- colnames(X)
  }
  
  K  <- length(Nks)
  P1 <- dim(X1)[2]
  P2 <- dim(X2)[2]
  
  # Unpack the constraints
  if( is.null(constraints$A) ){
    const <- create.constraints( P1, constraints)
    A     <- const$A
    Anull <- const$Anull
    B     <- const$B
  } else{
    A          <- constraints$A
    Anull      <- constraints$Anull
    B          <- constraints$B
    if( is.null(Anull) ){
      Anull <- create.constraints( P1, list(order="simple", node=1, decreasing=TRUE) )$Anull
    }
  }
  
  
  # Initialize values
  theta <- ginv( t(X)%*%X) %*% ( t(X)%*%Y )
  
  R   <- Y-X%*%theta
  ssq <- apply( as.matrix(1:K, nrow=1), 1 , 
                FUN=function(k, N1, N2, R ){ sum( R[N1[k]:N2[k]]^2 ) / Nks[k] } ,
                N1, N2, R)
  
  ssqvec  <- rep( ssq, Nks )  
    
  tsq  <- NULL
  
  theta1 <- theta
  ssq1   <- ssq
  tsq1   <- tsq  
  
  
  # Being the EM Algorithm convergence loop
  CONVERGE  <- 0
  iteration <- 0
  while( CONVERGE==0 ){
    iteration <- iteration+1
    
    R <- Y-X%*%theta
    
    if( verbose==TRUE ){
      message( "EM iteration " , iteration)
    }
      
    # Step 1: Estimate Sigma

      trace.vec <- (R/ssqvec)^2 - 1/ssqvec
    
    
    # trace.vec <- diag( PsiI%*%( R%*%t(R) )%*%PsiI - PsiI )
    ssq <- apply( as.matrix(1:K, nrow=1), 1 , 
                  FUN=function(k, ssq, Nks, N1, N2 , trv){
                    idx <- N1[k]:N2[k]
                    ssq[k] + ( (ssq[k]^2)/(Nks[k]) )*sum( trv[idx] )
                    } ,
                  ssq , Nks , N1 , N2 , trace.vec )
    
    ssqvec  <- rep( ssq, Nks)
        
    # Step 2a: Estimate Thetas
    # Update the blocks
    SiR   <- R / ssqvec      
      #theta[1:P1]  <- theta1[ 1:P1] + ginv(t(X1)%*%SigmaI%*%X1) %*% ((t(X1)%*%PsiI)%*%R )
      theta[1:P1]  <- theta1[1:P1] + ginv( t(X1)%*%(X1/ssqvec) ) %*% (t(X1) %*% SiR)
      
      
      if( is.null(X2)==FALSE ){
        #theta[(P1+1):(P1+P2)] <- ( theta1[ (P1+1):(P1+P2)] + 
        #                              ginv(t(X2)%*%SigmaI%*%X2) %*% ((t(X2)%*%PsiI)%*%R ) )
        X2SiR <- t(X2) %*% SiR
        theta[(P1+1):(P1+P2)] <- ( theta1[(P1+1):(P1+P2)] + 
                                     ginv(t(X2)%*%(X2/ssqvec)) %*% (X2SiR) )
      }
      
    
    cov.theta <- solve( t(X) %*% (X/ssqvec) )
    
    ## Apply order constraints / isotonization
    if( mySolver=="GLS"){
      wts <- solve( cov.theta )[1:P1, 1:P1, drop=FALSE]
    } else{
      wts <- diag( solve(cov.theta) )[1:P1]
    }
    theta[1:P1] <- activeSet(A, y = theta[1:P1], weights = wts, mySolver=mySolver )$x

    
    # Evaluate some convergence criterion
    rel.change <- abs(theta - theta1)/theta1
    if( mean(rel.change) < em.eps || iteration >= em.iter ){
      CONVERGE <- 1
    } else{
      theta1 <- theta
      ssq1   <- ssq
    }
    
  } # End converge loop (while)
  
  if( verbose==TRUE ){
    message("EM Algorithm ran for " , iteration , " iterations." )
  }
  
  theta        <- c(theta)
  names(theta) <- theta.names
    
  theta.null       <- theta
  theta.null[1:P1] <- activeSet( Anull, y = theta[1:P1], weights = wts , mysolver=mySolver )$x
  
  
  # Compute test statistic
  ts.glb <- tsf( theta=theta, theta.null=theta.null, cov.theta=cov.theta, B=B, A=A, Y=Y, X1=X1, 
                 X2=X2, U=U, tsq=tsq, ssq=ssq, Nks=Nks, Qs=Qs  )
  
  ts.ind <- tsf.ind(theta=theta, theta.null=theta.null, cov.theta=cov.theta, B=B, A=A, Y=Y, X1=X1, 
                    X2=X2, U=U, tsq=tsq, ssq=ssq, Nks=Nks, Qs=Qs )
  
  # Return the results
  em.results <- list(theta=theta, theta.null=theta.null, ssq=ssq, tsq=tsq,
                     cov.theta=cov.theta, ts.glb=ts.glb, ts.ind=ts.ind, mySolver=mySolver )
  
  return( em.results )

}

#' @description \code{clme_em_mixed} performs a constrained EM algorithm for linear mixed effects models.
#'
#' @rdname clme_em
#' @export
#'
clme_em_mixed <- function( Y, X1, X2 = NULL, U = NULL, Nks = dim(X1)[1],
                     Qs = dim(U)[2], constraints, mq.phi = NULL, tsf = lrt.stat, 
                     tsf.ind = w.stat.ind, mySolver="LS", em.iter = 500, 
                     em.eps =  0.0001, verbose = FALSE, ... ){
  
  
  if( verbose==TRUE ){
    message("Starting EM algorithm")    
  }
  
  N  <- sum(Nks)
  N1 <- 1 + cumsum(Nks) - Nks
  N2 <- cumsum(Nks)
  
  Q  <- length(Qs)
  Q1 <- 1 + cumsum(Qs) - Qs
  Q2 <- cumsum(Qs)
  
  X  <- as.matrix(cbind(X1, X2))
  theta.names <- NULL
  if( !is.null(colnames(X)) ){
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
  
  ssqvec  <- rep(  ssq,Nks)  
  
  
  if( is.null(mq.phi) ){
    mq.phi <- minque( Y=Y, X1=X1, X2=X2, U=U, Nks=Nks, Qs=Qs)
  }
  tsq    <- mq.phi[1:Q]
  tsqvec <- rep(tsq,Qs)
  
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
      SiR    <- R / ssqvec
      # U' * SigI * U
      U1         <- apply( U , 2 , FUN=function(x,sq){x*sq} , 1/sqrt(ssqvec) )
      tusu       <- t(U1) %*% U1
      diag(tusu) <- diag(tusu) + 1/tsqvec
      tusui      <- solve(tusu)
      PiR        <- SiR - (U %*% (tusui %*% (t(U)%*%SiR))) * (1/ssqvec)
      # PsiI
      BU    <- tusui %*% t(U)
      UBU   <- apply( as.matrix(1:sum(Nks)) , 1 , FUN=function(kk,uu,bu){ sum( uu[kk,]*bu[,kk] ) } , U, BU )
      SUBUS <- 1/ssqvec - UBU / ssqvec^2
      trace.vec <- PiR^2 - SUBUS      
    
    
    # trace.vec <- diag( PsiI%*%( R%*%t(R) )%*%PsiI - PsiI )
    ssq <- apply( as.matrix(1:K, nrow=1), 1 , 
                  FUN=function(k, ssq, Nks, N1, N2 , trv){
                    idx <- N1[k]:N2[k]
                    ssq[k] + ( (ssq[k]^2)/(Nks[k]) )*sum( trv[idx] )
                    } ,
                  ssq , Nks , N1 , N2 , trace.vec )
    
    ssqvec  <- rep(  ssq,Nks)
        
    # Step 2a: Estimate Thetas
    # Update the blocks
    SiR   <- R / ssqvec
    X1SiR <- t(X1) %*% SiR

    X1SiU <- t(X1) %*% (U/ssqvec)  #  X1SU
    USiR  <- t(U)  %*% SiR
    
    U1         <- apply( U , 2 , FUN=function(x,sq){x*sq} , 1/sqrt(ssqvec) )
    tusu       <- t(U1) %*% U1
    diag(tusu) <- diag(tusu) + 1/tsqvec
    tusui      <- solve(tusu)    
    
    
    #theta[1:P1]  <- theta1[ 1:P1] + ginv(t(X1)%*%SigmaI%*%X1) %*% ((t(X1)%*%PsiI)%*%R )
    theta[1:P1]  <- theta1[1:P1] + ginv( t(X1)%*%(X1/ssqvec) ) %*% (X1SiR - X1SiU%*%(tusui%*%USiR))
    
    
    if( is.null(X2)==FALSE ){
      #theta[(P1+1):(P1+P2)] <- ( theta1[ (P1+1):(P1+P2)] + 
      #                              ginv(t(X2)%*%SigmaI%*%X2) %*% ((t(X2)%*%PsiI)%*%R ) )
      X2SiU <- t(X2) %*% (U/ssqvec)
      X2SiR <- t(X2) %*% SiR
      theta[(P1+1):(P1+P2)] <- ( theta1[(P1+1):(P1+P2)] + 
                                   ginv(t(X2)%*%(X2/ssqvec)) %*% (X2SiR - X2SiU%*%(tusui%*%USiR)) )
    }
  
    
    
    
    # Step 2b: Estimate Tau
    # USiR  <- t(U)  %*% SiR # <-- previously calcualted
    USiU <- t(U)  %*% (U/ssqvec)
    UPiR <- USiR - USiU%*%(tusui%*%USiR)
    trace.vec.tau <- (UPiR^2) - diag(USiU) + diag( USiU%*%tusui%*%USiU )
    
    for( q in 1:Q ){
      # tau.idx <- Q1[q]:Q2[q]
      #Uq      <- as.matrix( U[,tau.idx] )
      #cq      <- dim(Uq)[2]
      #sumdg   <- sum(diag( t(Uq)%*%( PsiI%*%R%*%t(R)%*%PsiI - PsiI )%*%Uq ))      
      #tsq[q]  <- tsq1[q] + ((tsq1[q]^2)/cq)*sumdg
      tau.idx <- Q1[q]:Q2[q]
      cq      <- length(tau.idx)
      tsq[q]  <- tsq1[q] + ((tsq1[q]^2)/cq)*sum(trace.vec.tau[tau.idx])
    }
        
    XSiX      <- t(X) %*% (X/ssqvec)
    tsqvec    <- rep(tsq,Qs)  
    XSiU      <- t(X) %*% (U/ssqvec)
    cov.theta <- solve( XSiX - XSiU%*%tusui%*%t(XSiU) )
    
    ## Apply order constraints / isotonization
    if( mySolver=="GLS"){
      wts <- solve( cov.theta )[1:P1, 1:P1, drop=FALSE]
    } else{
      wts <- diag( solve(cov.theta) )[1:P1]
    }
    theta[1:P1] <- activeSet(A, y = theta[1:P1], weights = wts, mySolver=mySolver  )$x
    
    
    # Evaluate some convergence criterion
    rel.change <- abs(theta - theta1)/theta1
    if( mean(rel.change) < em.eps || iteration >= em.iter ){
      CONVERGE <- 1
    } else{
      theta1 <- theta
      ssq1   <- ssq
      tsq1   <- tsq
    }
    
  } # End converge loop (while)
  
  if( verbose==TRUE ){
    message("EM Algorithm ran for " , iteration , " iterations." )
  }
  
  theta        <- c(theta)
  names(theta) <- theta.names
  
  theta.null       <- theta
  theta.null[1:P1] <- activeSet( Anull, y = theta[1:P1], weights = wts , mySolver=mySolver )$x
  
  
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

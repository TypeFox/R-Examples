
.packageName <- "nets"

nets <- function( y , GN=TRUE , CN=TRUE , p=1 , lambda=stop("shrinkage parameter 'lambda' has not been set") , alpha.init=NULL , rho.init=NULL , algorithm='activeshooting' , weights='adaptive' , iter.in=100 , iter.out=2 , verbose=FALSE ){
  
  # input check
  if( !is.data.frame(y) & !is.matrix(y) ){
    stop("The 'y' parameter has to be a TxN matrix or a data.frame of data")
  }
  if( iter.in < 1 | iter.out < 1 ){
    stop("The 'iter' parameters have to be positive")
  }
  if( p<0 ){
    stop("The 'p' parameter has to be nonnegative")
  }
  if( !is.logical(verbose) ){
    stop("The 'verbose' parameter has to be TRUE or FALSE")
  }
  if( GN==FALSE & CN==FALSE ){
    stop("At least A or G have to be true")  
  }
  if( !(algorithm %in% c('shooting','activeshooting')) ){
    stop("The algorihm has to be set to: shooting' or 'activeshooting'")  
  }
  if( !(weights %in% c('adaptive','none')) ){
    stop("The lasso weights have to be set to: 'adaptive', or 'none'")  
  }
  
  # define variables
  T <- nrow(y)
  N <- ncol(y)
  P <- p
  if( !is.null(dimnames(y)[[2]]) ){
    labels <- dimnames(y)[[2]]
  } else {
    labels <- paste('V',1:N,sep='')
  }
  if( length(lambda)==1 ){
    lambda <- rep(lambda,2)
  }

  if( !is.null(alpha.init) ) alpha <- alpha.init
  else alpha <- rep(0,N*N*P)

  if( !is.null(rho.init) ) rho <- rho.init
  else rho <- rep(0,N*(N-1)/2)
  
  alpha.weights <- rep(1,N*N*P)
  rho.weights   <- rep(1,N*(N-1)/2)
  c.hat         <- 1/diag(cov(y))
  
  # ADAPTIVE WEIGHTS COMPUTATION
  if( weights=='adaptive' ){
    
    if( (T < N*P) ){ stop('Cannot compute adaptive weights by LS when T<N*P') }

    if( GN == TRUE ){
      
        x.aux <- matrix( 0 , T , N*P )
        for( p in 1:P ){
          x.aux[(p+1):T, ((p-1)*N+1):(p*N) ] <- y[1:(T-p),]
        }
        
        reg <- lm( y ~ 0+x.aux )
        A     <- coef(reg)
        eps   <- reg$resid
        
        alpha.pre <- c()
        for( p in 1:P ){
          alpha.pre <- c( alpha.pre , ( A[((p-1)*N+1):(p*N),] ) [1:(N*N)] )
        }

        alpha.weights <- 1/(abs(alpha.pre)+1e-4)

    } 
    else{
      eps      <- y
    }
    
    if( CN == TRUE ){
      C.hat       <- solve( cov(eps) )
      PC          <- -diag( c.hat**(-0.5) ) %*% C.hat %*% diag( c.hat**(-0.5) )	        
      rho.pre     <- PC[ upper.tri(PC) ]
      rho.weights <- 1/(abs(rho.pre)+1e-4)      
    }  
    else{
      iter.out <- 1
    }
  }
  
  if( algorithm == 'activeshooting' ){ algorithm <- 'nets_activeshooting' }
  else{ algorithm <- 'nets_shooting' }
  
  # call nets
  for( iter in 1:iter.out ){
    cat('iter',iter,'of',iter.out,'\n')
    run <- .C(algorithm,
             alpha        =as.double(alpha),
             rho          =as.double(rho), 
             alpha.weights=as.double(alpha.weights),
             rho.weights  =as.double(rho.weights),
             lambda       =as.double(lambda),
             y            =as.double(y),
             T            =as.integer(T),
             N            =as.integer(N),
             P            =as.integer(P),
             c.hat        =as.double(c.hat),
             GN           =as.integer(GN),
             CN           =as.integer(CN),
             v            =as.integer(verbose),
             m            =as.integer(iter.in),
	           rss          =as.double(0),
	           npar         =as.double(0))
  }
  
  # package results
  obj <- list()
  class(obj)    <- 'nets'
  obj$y         <- y
  obj$T         <- T
  obj$N         <- N
  obj$P         <- P
  obj$rss       <- run$rss/(N*T)
  obj$npar      <- run$npar
  obj$lambda    <- lambda
  obj$GN        <- GN
  obj$CN        <- CN
  
  if( GN == TRUE ){
    A.hat <- array(0,dim=c(N,N,P))
    for( p in 1:P ){
      for( i in 1:N ){
        A.hat[i,,p] <- run$alpha[ ((p-1)*N*N+(i-1)*N+1):((p-1)*N*N+i*N) ]
      }
    }
    dimnames(A.hat)[[1]] <- labels
    dimnames(A.hat)[[2]] <- labels    
    obj$A.hat     <- A.hat 
    obj$alpha.hat <- run$alpha 
  }
  else{
    obj$alpha.hat <- alpha
  }
  if( CN == TRUE ){
    C.hat <- matrix(0,N,N)
    for( i in 1:N ){
      C.hat[i,i] <- run$c[i]
      for( j in setdiff(1:N,i) ){
        c_ij       <- -run$rho[ (max(i,j)-1)*(max(i,j)-2)/2 + min(i,j) ] * sqrt( run$c.hat[i] * run$c.hat[j] )
        C.hat[i,j] <- c_ij
        C.hat[j,i] <- c_ij
      }
    }
    dimnames(C.hat)[[1]] <- labels
    dimnames(C.hat)[[2]] <- labels
    obj$C.hat     <- C.hat 
    obj$rho.hat   <- run$rho
    obj$c.hat     <- run$c.hat
  }
  else{
    obj$rho.hat   <- rho
    obj$c.hat     <- c.hat
  }
  
  # networks
  I.n <- diag(N)
  
  if( CN == TRUE ){
    PCN  <- C.hat > 0
    PCN[row(I.n) == col(I.n) ] <- 0
    obj$c.adj <- 1*(PCN > 0)
  }
  if( GN == TRUE ){
    DGN <- matrix(0,N,N)
    for( p in 1:P ){
      DGN <- DGN | A.hat[,,p]>0
    }
    DGN[row(I.n) == col(I.n) ] <- 0
    obj$g.adj    <- 1*DGN
  }
  if( CN == TRUE && GN==TRUE ){
    KL         <- t(I.n-DGN ) %*% C.hat %*% ( I.n-DGN )
    LRPCN      <- -diag( diag(KL)**(-0.5) ) %*% KL %*% diag( diag(KL)**(-0.5) )
    LRPCN[row(I.n) == col(I.n) ] <- 0
    obj$lr.adj <- 1*(LRPCN > 0)
  }
  
  return(obj)
}

print.nets <- function( x , ... ) {
   	cat( ' Time Series Panel Dimension: T=',x$T,' N=',x$N,'\n',sep='')
   	cat( ' VAR Lags P=',x$P,'\n',sep='')
    cat( ' RSS',x$rss,'Num Par',x$npar)
    cat( ' Lasso Penalty: ', x$lambda )
}

predict.nets <- function( object , newdata , ... ){

  # input check
  if( !is.data.frame(newdata) & !is.matrix(newdata) ){
    stop("The 'newdata' parameter has to be a TxN matrix or a data.frame of new observation")
  }
  
  x     <- object
  
  #
  T     <- nrow(newdata)
  N     <- x$N
  y.hat <- matrix(0,T,N)
  y     <- rbind(x$y[(x$T-x$P+1):x$T,],newdata)
  
  # call nets
  run <- .C( sprintf("nets_predict",algorithm),
             y.hat        =as.double(y.hat),
             y            =as.double(y),
             T            =as.integer(T),
             N            =as.integer(x$N),
             P            =as.integer(x$P),
             alpha        =as.double(x$alpha.hat),
             rho          =as.double(x$rho.hat),              
             c            =as.double(x$c.hat),
             GN           =as.integer(x$GN),
             CN           =as.integer(x$CN),
             rss          =as.double(0))
  
  y.hat <- matrix(run$y.hat,T,N)
  rss   <- run$rss/(T*N)
  
  # output
  list( y.hat=y.hat , rss=rss )
}
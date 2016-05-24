mlicaMAIN <-
function(prNCP, tol=0.0001, maxit=300, mu=1){# default values
  print("Entered MLica");
  #######################################################
  X <- prNCP$X;
  x <-  prNCP$x ;
  pEx <- prNCP$pEx ;
  pCorr <- prNCP$pCorr ;
  ntp <- dim(X)[1]; # number of "time points"
  ndim <- dim(X)[2]; # number of dimensions.
  ncp <- ncol(x); # number of components
  ######################################################
  # Estimated Source matrix
  Sest <- matrix(nrow=ntp,ncol=ncp);  
  # Choose initial random separating matrix B(=Ainv):
  B.old <- matrix(runif(ncp*ncp,0,1),nrow=ncp,ncol=ncp);
  B.o <- B.old ; # this is unwhitened B separating matrix needed to check convergence of algorithm
  # variables needed
  icount <- 0;  
  not.conv <- c(1,2); # arbitrary non-zero vector to give non-zero length
  y <- matrix(nrow=ntp,ncol=ncp);
  tmp <- matrix(nrow=ncp,ncol=ncp);
  beta <- vector(length=ncp);
  alpha <- vector(length=ncp);
  
  # STARTING ITERATIONS
  while ( (length(not.conv) > 0) && ( icount < maxit)  ){ # iteration loop : update B
   print(c("Entering iteration loop ",icount));
   # whiten y=B.old x first
   Cy <- B.old %*% t(B.old);
   svds <- eigen(Cy,symmetric=TRUE,EISPACK=FALSE) ;
   D <- diag(svds$values);
   E <- svds$vectors;
   Dinv <- solve(D, LINPACK=FALSE);
   # whitening matrix:
   V <- E %*% sqrt(Dinv) %*% t(E);

   # project B.old onto whitening set
   B.old <- V %*% B.old ;
  
   for ( g in 1:ntp){
    y[g,] <- B.old %*% x[g,]; # y is white
   }

   for ( c in 1:ncp){

     beta[c] <- 2*sum( y[,c]*tanh(y[,c]) )/ntp ;
     alpha[c] <- -1/(beta[c] -2 + 2* sum( tanh(y[,c])*tanh(y[,c]) )/ntp ) ;

     for ( c2 in 1:ncp){
       tmp[c,c2] <- -2*sum(tanh(y[,c])*y[,c2])/ntp ;
     }
       
   }
   print("Checkpt1");
   # Update separating matrix
   tmp <- diag(beta) + tmp ;
   B <- B.old + mu*diag(alpha) %*% tmp %*% B.old ;

   # Check convergence
   Dev <- abs(B - B.o);
   AvDev <- sum(Dev)/(ncp*ncp);
   print(c("AvDev=",AvDev));
         
   not.conv <- vector();
   not.conv <- as.vector( Dev[ Dev > tol ] );

   B.old <- B ;
   B.o <- B;
   icount <- icount + 1;

   for( g in 1:ntp ){
    Sest[g,] <- B %*% x[g,] ;
   }
   logL <- -2*sum( log( cosh(Sest) ) ) + ntp*log(abs(det(B)));

   print("iterated logL");
   print(logL);
   
 } # matches iteration loop  


  # find estimated source components
  # whiten y=Bx first (actually we only do PCA)(i.e variance is not scaled)
   Cy <- B %*% t(B);
   svds <- eigen(Cy,symmetric=TRUE,EISPACK=FALSE) ;
   D <- diag(svds$values);
   E <- svds$vectors;
   Dinv <- solve(D, LINPACK=FALSE);
   # whitening matrix (variance scaling here)
   V <- E %*% sqrt(Dinv) %*% t(E);

   # project B.old onto whitening set
   B <- V %*% B ;

  
  for( g in 1:ntp ){
    Sest[g,] <- B %*% x[g,] ;
  }
  
  Aest <- t(pEx %*% sqrt(pCorr) %*% t(B)) ; # this is now estimate of original mixing matrix, Ns x ncp matrix
  if ( length(not.conv) > 0 ){
    NotConv <- 1 ;
  }
  else { NotConv <- 0 ;}

  # logL computation
  logL <- -2*sum( log( cosh(Sest) ) ) + ntp*log(abs(det(B)));


  return(list(A=Aest,B=B,S=Sest,X=X,ncp=dim(Sest)[2],NC=NotConv,LL=logL));

}

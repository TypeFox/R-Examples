
######################
# Variational M step #
######################
VMstep<-function(X, nX, N, Q, tau, pCts0, directed=FALSE){


  pCts <- list()
  pCts$n <- pCts0$n + colSums(tau)
  
  pCts$eta <- t(tau)%*%X%*%tau
  if ( ! directed )
    pCts$eta[seq(1, Q^2, by=Q+1)] <- pCts$eta[seq(1, Q^2, by=Q+1)]/2

  pCts$eta <- pCts$eta + pCts0$eta
  
  pCts$zeta <- t(tau)%*%nX%*%tau
  if( ! directed )
    pCts$zeta[seq(1, Q^2, by=Q+1)] <- pCts$zeta[seq(1, Q^2, by=Q+1)]/2
    
  pCts$zeta <- pCts$zeta + pCts0$zeta
  
  return(pCts)
}             

###############################
# Variational E step - C call #
###############################
VEstep.C<-function( X, N, Q, mincut, maxcut, tau, pCts, dmaxInner, fpnbiter,
                    directed=FALSE) {
  
  dGamma <- 1
  nit    <- 1

  A <- digamma(pCts$eta)
  B <- digamma(pCts$zeta)
  C <- digamma(pCts$eta + pCts$zeta)
  
  if (directed)
    D <- B + t(B) - C - t(C)
  else
    D <- B - C
  
  E <- A - B

  digamma.n     <- digamma(pCts$n) 
  digamma.sum.n <- digamma(sum(pCts$n))

  CstMat <- matrix( rep(digamma.n - digamma.sum.n, N), N, Q, byrow=TRUE)

  i.directed <- directed
  
  res <- .C("VEstep_",
            as.integer(N),          # Number of nodes
            as.integer(Q),          # Number of classes
            as.double(X),           # Adjacency  matrix(N,N)
                                    # Remarks: <integer> converted to <double>
            as.double(D),           # Tau %*% D matrix(N,Q)
            as.double(E),           # Tau %*% E matrix(N,Q)
            as.double(CstMat),      # Constant matrix(N,Q)   
            as.double(mincut),      # Min. value cut-off   
            as.double(maxcut),      # Max. value cut-off   
            as.double(dmaxInner),   # Convergence criterion
            as.integer(fpnbiter),   # Max. iteration number
            as.integer(directed),   # Directed matrix
            res = as.double(tau),   # Taus (input/output)
            niter=as.integer(nit)   # Number of iterations
            )
  # cat(" [C, dir=", directed, " ] iter = ", res$niter, "\n")
  tau <- matrix(res$res, N, Q )
  return( tau )
}

###############
# Lower Bound #
###############
lowerBound<-function(tau, pCts0, pCts, directed=FALSE){

  Q <- dim(tau)[2]

  if ( directed ) {
    l <- - sum(tau*log(tau)) + lgamma(sum(pCts0$n)) - lgamma(sum(pCts$n)) - sum(lgamma(pCts0$n) - lgamma(pCts$n)) + sum(lbeta(pCts$eta, pCts$zeta) - lbeta(pCts0$eta, pCts0$zeta))
  } else {
    A <- lbeta(pCts$eta, pCts$zeta)
    B <- lbeta(pCts0$eta, pCts0$zeta)
    l <- - sum(tau*log(tau)) + lgamma(sum(pCts0$n)) - lgamma(sum(pCts$n)) - sum(lgamma(pCts0$n) - lgamma(pCts$n)) + (1/2)*sum(A[-seq(1, Q^2, by=Q+1)]) + sum(A[seq(1, Q^2, by=Q+1)]) - (1/2)*sum(B[-seq(1, Q^2, by=Q+1)]) - sum(B[seq(1, Q^2, by=Q+1)])
  }
  return(l)
}

#####################
# Variational Bayes #
#####################
VariationalBayes <-
  function(m, qmin, qmax, nbiter, fpnbiter, emeps, fpeps, directed=FALSE) {

  ## Edge matrix to Adjacency matrix
  N <- max(m)
  X<-matrix(0, N, N )
  X[ cbind( m[1,], m[2,]) ] <- 1
  if( ! directed )
    X[ cbind( m[2,], m[1,]) ] <- 1
  
  ## How to find the proper order.....
  ## readingOrder<-unique(as.numeric(m));

  
  if( is.null(qmax) )
    qmax <- qmin
  
  X <- sign(X) # check that the graph contains binary edges
  
  nX <- matrix(rep(1, N^2), N, N) - X
  nX[seq(1, N^2, by=N+1)] = 0
  

  vbOptions <- list()
  vbOptions$dmaxOuter <- emeps # Outer loop
  vbOptions$dmaxInner <- fpeps # Inner loop = Fixed point

  pCts0 <- list()

  # Parameters for ERMG Initilization 
  symetrize     <- ! directed
  loop          <- FALSE
  undirected    <- ! directed
  silent        <- TRUE
  nokmeans      <- TRUE
  kmeansnbclass <- 0
  kmeansnbiter  <- 30
  nbrEdges      <- length(m)/2         # size of the given array

  # Result index
  i.res <- 1                  
  y <-  vector("list", qmax-qmin+1)

  
  for (Q in qmin:qmax ) {

    maxcut <- log(.Machine$double.xmax) - log(Q)
    mincut <- log(.Machine$double.xmin)

    pCts0$n    <- rep(1, Q)
    pCts0$eta  <- matrix(rep(1, Q^2), Q, Q)
    pCts0$zeta <- matrix(rep(1, Q^2), Q, Q)

    # tau0 initialization
    
    # ERMG intialization parameters
    nbrNodes <- N
    nbrTaus  <- nbrNodes*Q

    xout <- .C("init_ermg",
          as.integer(symetrize),
          as.integer(loop),     
          as.integer(undirected),
          as.integer(silent),    
          as.integer(kmeansnbclass),
          as.integer(kmeansnbiter), 
          as.integer(nokmeans),     
          as.integer(Q),            
          as.integer(nbrEdges),
          as.integer(nbrNodes),
          as.integer(m),
          Taus      = double(nbrTaus)  
         )

    tau0 <- matrix( xout$Taus[1: (Q*nbrNodes)], Q,
                          nbrNodes,byrow=TRUE)

    # tau0[,readingOrder] <- tau0
    
    # Require a nbrNodes x Q matrix 
    tau <- t( tau0 )
  
    # Avoid overflow
    tau[tau<.Machine$double.xmin] <- .Machine$double.xmin

    delta <- 1
    i <- 0
    l <- list()


    while(  is.finite(delta) == TRUE
            & delta>vbOptions$dmaxOuter
            & (i < nbiter)
          ) {
      i <- i + 1
      pCts <- VMstep( X, nX, N, Q, tau, pCts0, directed )

      l[[i]] <- lowerBound(tau, pCts0, pCts, directed)
      tau <- VEstep.C( X, N, Q, mincut, maxcut, tau, pCts,
                       vbOptions$dmaxInner, fpnbiter, directed )
      
      if(i > 1) {
        delta <- abs((l[[i]] - l[[i-1]])/l[[i-1]])
      }
    }
  
    # Only intialization   
    if ( nbiter == 0 ) {
      pCts <- pCts0
      l[[1]] <- 0
    }
    
    # Store results
    # -------------
    #  
    # ICL item is nn fact the Bayesian criterion
    y[[i.res]]$criterion <- max( as.numeric(as.matrix(l)) ) 
    y[[i.res]]$alphas    <- pCts$n / sum( pCts$n ) 
    y[[i.res]]$Pis       <- pCts$eta / (pCts$eta+pCts$zeta)
    y[[i.res]]$Taus      <- t(tau)
    i.res <- i.res+1
  }
  
  return(y)  
   
}




  


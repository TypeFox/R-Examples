########################### Function to indentified parameters #################################################
identif <- function(A,
                     B,
                     P,
                     K=K,
                     trace.lba,
                     itmax.ide,
                     what){
 I <- nrow(A)
 J <- nrow(B)                    

 colnames(B) <- colnames(A) <- paste('LB',1:K,sep='')
 rownames(B) <- colnames(P)
 rownames(A) <- rownames(P)

 if(K==2){
    Pi <- A%*%t(B)

  #------------------------------------------------------------------------------
  #For K = 2 the algorithm described in van der Ark PhD thesis section 2.2 page 25
  #up to 29.
  #------------------------------------------------------------------------------
  #OUTER EXTREME BUDGETS
  if(what == 'outer'){

   mb <- apply(B, 
                1, 
                function(x) { 
                 t1 <- c(-x[2],x[1])/(x[1]-x[2])
                 Tb <- matrix(t1,ncol=2,byrow=T)
                 Tb%*%t(B)} )
   
  tiny <- 1e-14
  #mb is a matrix containing in each column the vector bj having a zero in 
  #its jth coordinate. Their end points S1j are in the line S1

  #The parametric equations of the line between S1 passing at end points of 
  #beta1 and beta2, having the following equations:
  #BETA1 + (BETA2 - BETA1)t, b[,1] + (b[,2]-b[,1])t. 
  #Now we find the value of t for which each bj #(coordinates in the axes) are 
  #greater than zero and then find which (b[,2]-b[,1]) are positive and negative.

  to <- -B[,1]/(B[,2]-B[,1])
  out1 <- to[to==max(to[to<0])]
  out2 <- to[to==min(to[to>0])]
  bout <- cbind(mb[,colnames(mb)==names(out1)], mb[,colnames(mb)==names(out2)])

  #The columns of bout are outer extreme budgets 
  #bt <- b[rownames(b)==colnames(bout),] 
  #Tb <- rbind(c(-bt[1,2],bt[1,1])/(bt[1,1]-bt[1,2]),
  #           c(-bt[2,2],bt[2,1])/(bt[2,1]-bt[2,2]))
  bi <- ginv(t(B))
  Tb <- t(bout)%*% bi
  #Tb is the transformation matrix to get the outer extreme budgets from b
  aout <- A%*%(solve(Tb))
  #Matrix of mixing parameters relative to outer extreme budgets
  
  colnames(aout) <- colnames(bout) <- colnames(B)

  rownames(aout) <- rownames(P)
  rownames(bout) <- colnames(P)

  iter_outer <- 'not applicable'

  res <- list(aout, 
              bout,
              iter_outer) 

  names(res) <- c("outer extreme mixing parameters" ,
                  "outer extreme budgets",
                  "identify iteration")
  return(res)      
  } else {
  #--------------------------------------------------------------------------
  #INNER EXTREME BUDGETS
  #the distances among all points representing the expected budgets are calculated
  # the two pints giving the max of those distances are the inner extreme budgets
  di <- ni <- vector()
  v <- 1
  for(i in 1:(I-1)){ for(j in (i+1):I) {
   ni[v] <- paste(i,".",j, sep="")
   di[v] <- sum((Pi[i,]-Pi[j,])^2)
   v <- v + 1 
  }  }
  names(di) <- ni 
  na <- names(which(di==max(di)))[1]
  pl <- as.numeric(unlist(strsplit(na, "\\.")))
  binn <- t(Pi[pl,]) #The rows of binn are the inner extreme budgets
  names(dimnames(binn)) <- NULL
  #Ta <- a[rownames(a)==colnames(binn),]
  bi <- ginv(t(B))
  Ta <- t(binn)%*% bi
  #Ta is the transformation matrix to get the inner extreme budgets from b
  ainn <- A%*%(solve(Ta))
  #Matrix of mixing parameters relative to inner extreme budgets
  #---------------------------------------------------------------------------

  colnames(ainn)  <- colnames(binn) <- colnames(B)

  rownames(ainn) <- rownames(P)
  rownames(binn) <- colnames(P)

  iter_inner <- 'not applicable'

  res <- list(ainn, 
              binn,
              iter_inner)

  names(res) <- c("inner extreme mixing parameters" ,
                  "inner extreme budgets",
                  "identify iteraction")
  return(res)   
  }
 } else {

  Q <- choose(K, 2)
  Qc <- combn(K, 2)
  P[P==0] <- 1e-6
  B[B==0] <- 1e-9
  A[A==0] <- 1e-9
  tn0 <- as.vector(diag(K))
  # The function to be minimized for identification

  if(what == 'outer'){
  chisq <- function(tn, P, A, B, K, Qc, Q, I, J){
   T1 <- matrix(tn, K, byrow=T)
   Bn1 <- T1 %*% t(B)
   Dchi <- rep(0, Q)
   for (q in 1:Q) {
    for (j in 1:J) Dchi[q] <- ((Bn1[Qc[1, q], j] - Bn1[Qc[2, q], j])^2/sum(P[, j]))  + Dchi[q]
   }
   sDchi <- 1/sum(sqrt(Dchi)) # outer budgets
  }

  } else {
  chisqI <- function(tn, P, A, B, K, Qc, Q, I, J){
   T1 <- matrix(tn, K, byrow=T)
   Bn1 <- T1 %*% t(B)
   Dchi <- rep(0, Q)
   for (q in 1:Q) {
    for (j in 1:J) Dchi[q] <- ((Bn1[Qc[1, q], j] - Bn1[Qc[2, q], j])^2/sum(P[, j]))  + Dchi[q]
   }
   sDchi <- sum(sqrt(Dchi))   # inner budgets
  }
  }
  # The equallity constrained function , on matrix T.
  heq <- function(tn, P, A, B, K, Qc, Q, I, J)  {
         h <- rowSums(matrix(tn, K, byrow=T)) - rep(1, K)
  }

  # The inequallity constrained function , on matrix T.
   hin <- function(tn, P, A, B, K, Qc, Q, I, J)  {
         T1 <- matrix(tn, K, byrow=T)
         Bn1 <- T1 %*% t(B)
         An1 <- A %*% ginv(T1)
         VB <- c(as.vector(Bn1), as.vector(An1))
         h <- rep(0, 1)
         l <- J*K + I*K
         for(i in 1 : l) h[i] <- VB[i] + 0.001
        h
        }
  
  if(what == 'outer'){
  itmax.ala <- round(0.1*itmax.ide)
  itmax.opt <- round(0.9*itmax.ide)
  # Identification of the latent parameters
  Topt <- constrOptim.nl(par = tn0, 
                         fn = chisq,
                         P=P,
                         A=A,
                         B=B,
                         K=K,
                         Qc=Qc,
                         Q=Q,
                         I=I,
                         J=J,
                         heq = heq, 
                         hin= hin,
                         control.outer = list(trace=trace.lba,
                                              itmax=itmax.ala),
                         control.optim = list(maxit=itmax.opt))

  bout <- matrix(Topt$par, ncol=K, byrow=TRUE) %*% t(B)
  aout <- t(A %*% ginv(matrix(Topt$par, ncol=K, byrow=TRUE)))

  bout[which(bout < 0)] <- 0
  aout[which(aout < 0)] <- 0

  rownames(aout) <- rownames(bout) <- colnames(A)

  iter_ide <- round(as.numeric(Topt$counts[2]/Topt$outer.iterations))+Topt$outer.iterations

  res <- list(t(aout), 
              t(bout),
              iter_ide) 

  names(res) <- c("outer extreme mixing parameters" ,
                  "outer extreme budgets",
                  "identify iteraction")
  return(res)    

  } else {

  itmax.alai <- round(.1*itmax.ide)
  itmax.opti <- round(.9*itmax.ide)

  ToptI <- constrOptim.nl(par = tn0, 
                          fn = chisqI,
                          P=P,
                          A=A,
                          B=B,
                          K=K,
                          Qc=Qc,
                          Q=Q,
                          I=I,
                          J=J,
                          heq = heq, 
                          hin= hin,
                          control.outer = list(trace=trace.lba,
                                               itmax=itmax.alai),
                          control.optim = list(maxit=itmax.opti))

  binn <- matrix(ToptI$par, ncol=K, byrow=TRUE) %*% t(B)
  ainn <- t(A %*% ginv(matrix(ToptI$par, ncol=K, byrow=TRUE)))

  binn[which(binn < 0)] <- 0
  ainn[which(ainn < 0)] <- 0
  
  rownames(ainn) <- rownames(binn) <- colnames(A)

  iter_ide <- round(as.numeric(ToptI$counts[2]/ToptI$outer.iterations))+ToptI$outer.iterations

  res <- list(t(ainn), 
              t(binn),
              iter_ide)

  names(res) <- c("inner extreme mixing parameters" ,
                  "inner extreme budgets",
                  "identify iteration")
  return(res)    
  
  }
 }
}

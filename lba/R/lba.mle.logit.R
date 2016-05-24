lba.mle.logit <- function(obj        ,
                          A          ,
                          B          ,
                          K          ,
                          cA         ,
                          cB         ,
                          logitA     , 
                          logitB     ,  
                          omsk       ,
                          psitk      ,
                          S          ,
                          T          ,
                          itmax.ide  ,
                          trace.lba  ,
                          ...)
{


 #===========================================================================
 #Multinomial logit on both mixing parameters and latent components
 #===========================================================================
 if(all(c(!is.null(logitA), !is.null(logitB)))){

  results <- lba.mle.logit.AB(obj         = obj,
                              logitA      = logitA, 
                              logitB      = logitB,
                              omsk        = omsk,
                              psitk       = psitk,
                              K           = K,
                              S           = S,
                              T           = T,
                              itmax.ide   = itmax.ide,
                              trace.lba   = trace.lba,
                              ...) 
 } else 

  #===========================================================================
  #Multinomial logit on mixing parameters only
  #===========================================================================
  if(all(c(!is.null(logitA), is.null(logitB)))){

   results <- lba.mle.logit.A(obj         = obj,
                              B           = B,
                              cB          = cB,
                              logitA      = logitA, #design matrix for row covariates IxS
                              omsk        = omsk,
                              K           = K,
                              S           = S,
                              itmax.ide   = itmax.ide,
                              trace.lba   = trace.lba,
                              ...)
  } else {

   #===========================================================================
   #Multinomial logit on latent components only
   #===========================================================================
   results <- lba.mle.logit.B(obj       = obj,
                              A         = A,
                              cA        = cA,
                              logitB    = logitB, #design matrix for row covariates IxS 
                              psitk     = psitk,
                              K         = K,                            
                              T         = T,
                              itmax.ide = itmax.ide,
                              trace.lba = trace.lba,
                              ...)
  }
}

######################## AUXILIAR FUNCTIONS ##############################
lba.mle.logit.AB <- function(obj      ,
                             logitA   , 
                             logitB   ,
                             omsk     ,
                             psitk    ,
                             K        ,
                             S        ,
                             T        ,
                             itmax.ide,
                             trace.lba,
                             ...)
{

 I  <- nrow(obj)       # row numbers of data matrix
 J  <- ncol(obj)       # column numbers of data matrix

 #-----------------------------------------------------------------------------
 if (is.null(omsk)) {
  omsk <- matrix(c(rep(0,S), rnorm(S*(K-1))), ncol = K) }

 if (is.null(psitk)) {
  psitk <- matrix(rnorm(T*K), ncol = K) }

 #----------------------------------------------------------------------------
 #third case multinomial logit constraints on mixing parameters, and
 #on the latent budgets.
 #---------------------------------------------------------------------------

 mw <- function(xx, 
                obj,
                K, 
                I, 
                J,
                logitA,
                logitB,
                S,
                T){

  omsk <- matrix(xx[(1):(S*K)], ncol = K)
  psitk <- matrix(xx[(S*K+1):(S*K+T*K)], ncol = K)

  # creating A from omsk (om(s,k) )
  A <- matrix(0,nrow=I,ncol=K)
  for(i in 1:I){
   for(k in 1:K){
    for(n in 1:K){
     a <- 1
     for(s in 1:S){
      if(exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))==Inf){
       a <- a*1e6
     }else{ a <- a*exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))}  }
     A[i,k] <- A[i,k] + a   }
    A[i,k] <- 1/A[i,k] }  }

  A <- A/rowSums(A)

  # A is a IxK matrix

  #creating  B from psitk
  B <- matrix(0,nrow=J,ncol=K)
  for(j in 1:J){
   for(k in 1:K){
    for(n in 1:J){
     a <- 1
     for(t in 1:T){
      if(exp(psitk[t,k]*(logitB[n,t]-logitB[j,t]))==Inf) {
       a <- a*1e6
     }else{ a <- a*exp(psitk[t,k]*(logitB[n,t]-logitB[j,t])) }  }
     B[j,k] <- B[j,k] + a  }
    B[j,k] <- 1/B[j,k] }  }

  B <- t(t(B)/colSums(B))
  # B is a JxK matrix

  pjki <- rep(0,I*J*K) #this will become pijk/pi+ see van der Ark page 80
  m <- 0
  # this makes i=1,j=1,k=1; i=2,j=1,k=1; i=3,j=1,k=1...; i=I,j=1,k=1 and so on.
  for(k in 1:K) for(j in 1:J) for(i in 1:I) { m <-  m +1
  pjki[m] <- A[i,k]*B[j,k] }

  pip <- rowSums(obj)/sum(obj)
  pjki[pjki <=0] <- 1e-7
  mi <- matrix(pjki, nrow=I)
  pijk <- mi*pip #this is pijk see van der Ark page 80
  nijk <- as.vector(pijk*sum(obj))
  mw <- -sum(nijk * log(pjki))  #-loglikelihood function 

 }

 x0 <- c(as.vector(omsk), as.vector(psitk))

 # finding the ls estimates

 xab <- optim(par    = x0,
              fn      = mw,
              obj     = obj,
              K       = K,
              I       = I,
              J       = J,
              logitA  = logitA,
              logitB  = logitB,
              S       = S,
              T       = T,
              method  = 'BFGS',
              control = list(trace=trace.lba,
                             maxit=itmax.ide))

 omsk  <- matrix(xab$par[1:(S*K)], ncol = K)
 psitk <- matrix(xab$par[(S*K+1):(S*K+T*K)], ncol = K)

 #creating  A from omsk
 A <- matrix(0,nrow=I,ncol=K)
 for(i in 1:I){
  for(k in 1:K){
   for(n in 1:K){
    a <- 1
    for(s in 1:S){
     if(exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))==Inf){
      a <- a*1e6
    }else{ a <- a*exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))}  }
    A[i,k] <- A[i,k] + a   }
   A[i,k] <- 1/A[i,k] }  }

 A <- A/rowSums(A)

 #creating  B from psitk
 B <- matrix(0,nrow=J,ncol=K)
 for(j in 1:J){
  for(k in 1:K){
   for(n in 1:J){
    a <- 1
    for(t in 1:T){
     if(exp(psitk[t,k]*(logitB[n,t]-logitB[j,t]))==Inf) {
      a <- a*1e6
    }else{ a <- a*exp(psitk[t,k]*(logitB[n,t]-logitB[j,t])) }  }
    B[j,k] <- B[j,k] + a  }
   B[j,k] <- 1/B[j,k] }  }

 B <- t(t(B)/colSums(B))
 # B is a JxK matrix

 B <- t(t(B)/colSums(B))
 # B is a JxK matrix

 colnames(A) <- colnames(B) <- colnames(omsk) <- colnames(psitk) <- paste('LB',1:K,sep='')

 pimais <- rowSums(obj)/sum(obj)

 pk <- pimais %*% A # budget proportions

 P <- obj/rowSums(obj)  
 rownames(A) <- rownames(P)
 rownames(B) <- colnames(P)

 pij <- A %*% t(B) # expected budget 

 residual <- P - pij

 val_func <- xab$value

 iter_ide <- as.numeric(xab$counts[2])

 rescB <- rescaleB(obj,
                   A,
                   B)

 colnames(rescB) <- colnames(B)
 rownames(rescB) <- rownames(B)

 results <- list(P,
                 pij,
                 residual,
                 A, 
                 B,
                 rescB,
                 pk,
                 val_func,
                 iter_ide,
                 omsk,
                 psitk)

 names(results) <- c('P',
                     'pij',
                     'residual',
                     'A',
                     'B',
                     'rescB',
                     'pk',
                     'val_func',
                     'iter_ide',
                     'omsk',
                     'psitk')

 #  names(results) <- c('Composition data matrix',
 #                      'Expected budget',
 #                      'Residual matrix',
 #                      'Mixing parameters',
 #                      'Latent budgets',
 #                      'Budget proportions',
 #                      'Value of the -loglik function',
 #                      'Number of iteractions',
 #                      'Logit parameters omsk',
 #                      'Logit parameters psitk')
 # 
 class(results)  <- c("lba.mle.logit",
                      "lba.mle")
 invisible(results) 

} 

lba.mle.logit.A <- function(obj       ,
                            B        ,
                            cB       ,
                            logitA   ,
                            omsk     ,
                            K        ,
                            S        ,
                            itmax.ide,
                            trace.lba,
                            ...)

{
 #The matrices caki and cbjk contain the constraint values of the mixing
 #parameters and latent components respectively.
 #For fixed value constraint use the values at respective location in the matrix.
 #For aki, all row sums must be less or equal 1. For bjk all column sums must be
 #less or equal 1.
 #For equality value constraint use whole numbers starting form 2. Same numbers
 #at diffferent locations of the matrix show equal parameters.
 #USE NA TO FILL UP THE REST OF THE MATRICES.

 I  <- nrow(obj)       # row numbers of data matrix
 J  <- ncol(obj)       # column numbers of data matrix

 #-----------------------------------------------------------------------------
 if (is.null(omsk)) {
  omsk <- matrix(c(rep(0,S), rnorm(S*(K-1))), ncol = K) }

 #BUILDING B
 if(!is.null(cB) & is.null(B)){
  B <- t(constrainAB(t(cB)))
  } else  { if(is.null(cB) & is.null(B)){
   #creating random generated values for beta(j|k)
   B <- t(rdirich(K, runif(J))) }  }

 #============================================================================

 #============================================================================

 #----------------------------------------------------------------------------
 #first case multinomial logit constraints on mixing parameters, but not
 #in the latent budgets.
 #---------------------------------------------------------------------------

 mw <- function(xx, 
                obj,
                cB,
                K, 
                I, 
                J,
                logitA,
                S){

  y <- length(xx)- (J*K+S*K)
  omsk <- matrix(xx[(y+J*K +1):(y+J*K+S*K)], ncol = K)
  B <- matrix(xx[(y+1):(y+J*K)], ncol = K)

  #creating A from omsk
  A <- matrix(0,nrow=I,ncol=K)
  for(i in 1:I){
   for(k in 1:K){
    for(n in 1:K){
     a <- 1
     for(s in 1:S){
      if(exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))==Inf){
       a <- a*1e6
     }else{ a <- a*exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))}  }
     A[i,k] <- A[i,k] + a   }
    A[i,k] <- 1/A[i,k] }  }

  A <- A/rowSums(A)

  # A is a IxK matrix

  if(!is.null(cB)) {
   mincb <- min(cB, na.rm=TRUE)
   if(mincb < 1){
    posFb   <- which(cB<1, arr.ind = T)
    B[posFb] <- cB[posFb]  }   }

  pjki <- rep(0,I*J*K) #this will become pijk/pi+ see van der Ark page 80
  m <- 0
  # this makes i=1,j=1,k=1; i=2,j=1,k=1; i=3,j=1,k=1...; i=I,j=1,k=1 and so on.
  for(k in 1:K) for(j in 1:J) for(i in 1:I) { m <-  m +1
  pjki[m] <- A[i,k]*B[j,k] }

  pip <- rowSums(obj)/sum(obj)                                          
  pjki[pjki <=0] <- 1e-7
  mi <- matrix(pjki, nrow=I)
  pijk <- mi*pip #this is pijk see van der Ark page 80
  nijk <- as.vector(pijk*sum(obj))
  mw <- -sum(nijk * log(pjki))  #-loglikelihood function 
 }

 #============================================================================
 #         heq function
 #============================================================================
 heq <- function(xx, 
                 obj,
                 cB,
                 K, 
                 I, 
                 J,
                 logitA,
                 S){

  # construction of matrices omsk and B(B) from input vector x
  y <- length(xx)- (J * K + S * K )
  omsk <- matrix(xx[(y+J*K +1):(y+J*K+S*K)], ncol = K)
  B <- matrix(xx[(y+1):(y+J*K)], ncol = K)

  # creating random generated values from om(s,k)
  A <- matrix(0,nrow=I,ncol=K)

  for(i in 1:I){

   for(k in 1:K){

    for(n in 1:K){

     a <- 1

     for(s in 1:S){

      if(exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))==Inf){

       a <- a*1e6

     }else{ a <- a*exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))}  }

     A[i,k] <- A[i,k] + a   }

    A[i,k] <- 1/A[i,k] }  }

  A <- A/rowSums(A)


  if(!is.null(cB)) {
   mincb <- min(cB, na.rm=TRUE)
   maxcb <- max(cB, na.rm=TRUE) 
   if(maxcb > 1){
    bl <- list()
    for(i in 2:max(cB, na.rm=TRUE)){
     bl[[i-1]] <- which(cB==i, arr.ind=TRUE)
     bl[[i-1]] <- bl[[i-1]][order(bl[[i-1]][,1],bl[[i-1]][,2]),]
    }  

    #this code forces the corresponding betasjk's to be equal
    h <- 0
    b <- 0
    e <- 0
    for(j in 1:(maxcb-1)) {
     b[j] <- nrow(bl[[j]])
     for(i in 1:(b[j])){  e <- e+1
     h[e]<- B[bl[[j]][i,1],bl[[j]][i,2]]- xx[j] 
   }
   }  
   } 

   #this code forces the fixed constraints to be preserved.
   if(mincb < 1) {
  
    posFb   <- which(cB<1, arr.ind = T)
    
    if(maxcb > 1){
    
     h[(length(h)+1):(length(h)+ nrow(posFb))] <-(B[posFb] - cB[posFb]) 
    
    }else{
    
     h <- 0  
     h[1:nrow(posFb)] <- (B[posFb] - cB[posFb]) 
  } 
  } 
  }

  if(is.null(cB)){
  
   h <- c(colSums(B) - rep(1,(K)), xx[(y+J*K +1):(y+J*K+S)])
  
  }else{
  
   h[(length(h)+1):(length(h)+ K + S)] <- c((colSums(B) - rep(1,(K))), xx[(y+J*K +1):(y+J*K+S)])  
  }   
  h
 }

 #=========================================================================
 #                hin  function
 #=========================================================================
 hin <- function(xx, 
                 obj,
                 cB,
                 K, 
                 I, 
                 J,
                 logitA,
                 S){
  y <- length(xx)- (J*K + S*K )
  h <- xx[(y+1):(y+J*K)] + 1e-7
  h
 }

 #===========================================================================
 #===========================================================================

 if(!is.null(cB)){
  maxcb <- max(cB, na.rm=TRUE) 
  if(maxcb > 1) { #there are equality parameters in cB
   #list containing in each element the positions of the equality parameters
   #of matrix B(B) that are equal among them.
   bl <- list()
   for(i in 2:max(cB, na.rm=TRUE)){
    bl[[i-1]] <- which(cB==i, arr.ind=TRUE)
    bl[[i-1]] <- bl[[i-1]][
                           order(bl[[i-1]][,1],bl[[i-1]][,2]),]}
   m <- sum(sapply(bl, function(x) 1))
   a <- rep(0,m)
   } else { m <- 0 
   a <- rep(0,m) }
  } else { m <- 0 
  a <- rep(0,m) }

 x0 <- c(a, as.vector(B), as.vector(omsk))

 # finding the mle estimates

 itmax.ala <- round(.1*itmax.ide)
 itmax.opt <- round(.9*itmax.ide)

 xab <- constrOptim.nl(par         = x0,
                       fn          = mw,
                       cB          = cB,
                       logitA      = logitA,
                       obj          = obj,
                       K           = K,
                       I           = I,
                       J           = J,
                       S           = S,
                       heq         = heq,
                       hin         = hin,
                       control.outer=list(trace=trace.lba,
                                          itmax=itmax.ala),
                       control.optim=list(maxit=itmax.opt))

 y <- length(xab$par)- (J * K + S * K )
 omsk <- matrix(xab$par[(y+J*K +1):(y+J*K+S*K)], ncol = K)
 B <- matrix(xab$par[(y+1):(y+J*K)], ncol = K)

 A <- matrix(0,nrow=I,ncol=K)
 for(i in 1:I){
  for(k in 1:K){
   for(n in 1:K){
    a <- 1
    for(s in 1:S){
     if(exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))==Inf){
      a <- a*1e6
    }else{ a <- a*exp(logitA[i,s]*(omsk[s,n]-omsk[s,k]))}  }
    A[i,k] <- A[i,k] + a   }
   A[i,k] <- 1/A[i,k] }  }

 pimais <- rowSums(obj)/sum(obj)

 P <- obj/rowSums(obj)   
 pk <- pimais %*% A # budget proportions 

 colnames(A) <- colnames(B) <- colnames(omsk) <- colnames(pk) <- paste('LB',1:K,sep='') 

 pij <- A %*% t(B) # expected budget

 rownames(A) <- rownames(P)
 rownames(B) <- colnames(P)
 
 residual <- P - pij

 val_func <- xab$value

 iter_ide <- round(as.numeric(xab$counts[2]/xab$outer.iterations)) + xab$outer.iterations

 rescB <- rescaleB(obj,
                   A,
                   B)

 colnames(rescB) <- colnames(B)
 rownames(rescB) <- rownames(B)

 results <- list(P,
                 pij,
                 residual,
                 A, 
                 B,
                 rescB,
                 pk,
                 val_func,
                 iter_ide,
                 omsk)

 names(results) <- c('P',
                     'pij',
                     'residual',
                     'A',
                     'B',
                     'rescB',
                     'pk',
                     'val_func',
                     'iter_ide',
                     'omsk')

 #  names(results) <- c('Composition data matrix',
 #                      'Expected budget',
 #                      'Residual matrix',
 #                      'Mixing parameters',
 #                      'Latent budgets',
 #                      'Budget proportions',
 #                      'Value of the ls function',
 #                      'Number of iteractions',
 #                      'Logit parameters omsk')
 # 
 class(results)  <- c("lba.mle.logit",
                      "lba.mle")
 invisible(results)

} 

lba.mle.logit.B <- function(obj       ,
                            A        ,
                            cA       ,
                            logitB   , 
                            psitk    ,
                            K        ,  
                            T        ,
                            itmax.ide,
                            trace.lba,
                            ...)
{
 #The matrices caki and cbjk contain the constraint values of the mixing
 #parameters and latent components respectively.
 #For fixed value constraint use the values at respective location in the matrix.
 #For aki, all row sums must be less or equal 1. For bjk all column sums must be
 #less or equal 1.
 #For equality value constraint use whole numbers starting form 2. Same numbers
 #at diffferent locations of the matrix show equal parameters.
 #USE NA TO FILL UP THE REST OF THE MATRICES.
 I  <- nrow(obj)       # row numbers of data matrix
 J  <- ncol(obj)       # column numbers of data matrix

 #-----------------------------------------------------------------------------
 if (is.null(psitk)) {
  psitk <- matrix(rnorm(T*K), ncol = K) }

 #BUILDING cA
 if(!is.null(cA) & is.null(A)){
  A <- constrainAB(cA)
  } else  { if(is.null(cA) & is.null(A)){
   #creating random generated values for beta(j|k)
   A <- rdirich(I, runif(K)) }  }

 #============================================================================

 #----------------------------------------------------------------------------
 #second case multinomial logit constraints on latent budgets, but not
 #in the mixing parameters.
 #---------------------------------------------------------------------------

 mw <- function(xx, 
                obj,
                cA,
                K, 
                I, 
                J,
                logitB,
                T){

  y <- length(xx)- (T*K+I*K)
  A <- matrix(xx[(y+1):(y+I*K)], ncol = K)
  psitk <- matrix(xx[(y+I*K+1):(y+I*K+T*K)], ncol = K)

  #creating  B from psitk
  B <- matrix(0,nrow=J,ncol=K)
  for(j in 1:J){
   for(k in 1:K){
    for(n in 1:J){
     a <- 1
     for(t in 1:T){
      if(exp(psitk[t,k]*(logitB[n,t]-logitB[j,t]))==Inf) {
       a <- a*1e6
     }else{ a <- a*exp(psitk[t,k]*(logitB[n,t]-logitB[j,t])) }  }
     B[j,k] <- B[j,k] + a  }
    B[j,k] <- 1/B[j,k] }  }

  B <- t(t(B)/colSums(B))
  # B is a JxK matrix


  if(!is.null(cA)) {
   minca <- min(cA, na.rm=TRUE)
   if(minca < 1){
    posFa   <- which(cA<1, arr.ind = T)
    A[posFa] <- cA[posFa]  }   }

  pjki <- rep(0,I*J*K) #this will become pijk/pi+ see van der Ark page 80
  m <- 0
  # this makes i=1,j=1,k=1; i=2,j=1,k=1; i=3,j=1,k=1...; i=I,j=1,k=1 and so on.
  for(k in 1:K) for(j in 1:J) for(i in 1:I) { m <-  m +1
  pjki[m] <- A[i,k]*B[j,k] }

  pip <- rowSums(obj)/sum(obj)
  pjki[pjki <=0] <- 1e-7
  mi <- matrix(pjki, nrow=I)
  pijk <- mi*pip #this is pijk see van der Ark page 80
  nijk <- as.vector(pijk*sum(obj))
  mw <- -sum(nijk * log(pjki))  #-loglikelihood function 
 }

 #============================================================================
 #         heq function
 #============================================================================
 heq <- function(xx,
                 obj,
                 cA,
                 K, 
                 I, 
                 J,
                 logitB,
                 T){

  # construction of matrices A(A) and B(B) from input vector x

  y <- length(xx)- (T*K+I*K)
  A <- matrix(xx[(y+1):(y+I*K)], ncol = K)
  psitk <- matrix(xx[(y+I*K+1):(y+I*K+T*K)], ncol = K)

  B <- matrix(0,nrow=J,ncol=K)
  for(j in 1:J){
   for(k in 1:K){
    for(n in 1:J){
     a <- 1
     for(t in 1:T){
      if(exp(psitk[t,k]*(logitB[n,t]-logitB[j,t]))==Inf) {
       a <- a*1e6
     }else{ a <- a*exp(psitk[t,k]*(logitB[n,t]-logitB[j,t])) }  }
     B[j,k] <- B[j,k] + a  }
    B[j,k] <- 1/B[j,k] }  }

  B <- t(t(B)/colSums(B))
  # B is a JxK matrix


  if(!is.null(cA)) {
   minca <- min(cA, na.rm=TRUE)
   maxca <- max(cA, na.rm=TRUE)
   if(maxca > 1){
    al <- list()
    for(i in 2:max(cA, na.rm=TRUE)){
     al[[i-1]] <- which(cA==i, arr.ind=TRUE)
     al[[i-1]] <- al[[i-1]][
                            order(al[[i-1]][,1],al[[i-1]][,2]),] }

    #this code forces the corresponding alphasik's to be equal
    h <- 0
    a <- 0
    e <- 0
    for(j in 1:(maxca-1)) {
     a[j] <- nrow(al[[j]])
     for(i in 1:(a[j])){ e <- e+1
     h[e]<-
      A[al[[j]][i,1],al[[j]][i,2]]- xx[j] } } } 

   #this code forces the fixed constraints to be preserved.
   if(minca < 1) {
    posFa   <- which(cA<1, arr.ind = T)
    if(maxca > 1){
     h[(length(h)+1):(length(h)+ nrow(posFa))] <- 
      (A[posFa] - cA[posFa]) 
    }else{
     h <- 0
     h[1:nrow(posFa)] <- (A[posFa] - cA[posFa]) }  } }

  if(is.null(cA)){
   h <- c(rowSums(A) - rep(1,I))
  }else{
   h[(length(h)+1):(length(h)+ I)] <- c(rowSums(A) - rep(1,I))  }
  h
 }

 #=========================================================================
 #                hin  function
 #=========================================================================
 hin <- function(xx, 
                 obj,
                 cA,
                 K, 
                 I, 
                 J,
                 logitB,
                 T){
  y <- length(xx)- (I*K + T*K )
  h <- xx[(y+1):(y+I*K)] + 1e-7
  h
 }

 #===========================================================================
 #===========================================================================
 if(!is.null(cA)) { #there are equality parameters in cA
  #list containing in each element the positions of the equality parameters
  #of matrix A (A) that are equal among them.
  maxca <- max(cA, na.rm=TRUE)
  if(maxca > 1){
   al <- list()
   for(i in 2:max(cA, na.rm=TRUE)){
    al[[i-1]] <- which(cA==i, arr.ind=TRUE)
    al[[i-1]] <- al[[i-1]][order(al[[i-1]][,1],al[[i-1]][,2]),] }
   m <- sum(sapply(al, function(x) 1))
   a <- rep(0,m)
   } else { m <- 0
   a <- rep(0,m) }
  } else { m <- 0
  a <- rep(0,m) }

 x0 <- c(a, as.vector(A), as.vector(psitk))

 # finding the ls estimates

 itmax.ala <- round(.1*itmax.ide)
 itmax.opt <- round(.9*itmax.ide)

 xab <- constrOptim.nl(par     = x0,
                       fn      = mw,
                       cA      = cA,
                       logitB  = logitB,
                       obj      = obj,
                       K       = K,
                       I       = I,
                       J       = J,
                       T       = T,
                       heq     = heq,
                       hin     = hin,
                       control.outer=list(trace=trace.lba,
                                          itmax=itmax.ala),
                       control.optim=list(maxit=itmax.opt))
 # 
 y <- length(xab$par)- (T*K+I*K)
 psitk <- matrix(xab$par[(y+I*K+1):(y+I*K+T*K)], ncol = K)
 A <- matrix(xab$par[(y+1):(y+I*K)], ncol = K)

 B <- matrix(0,nrow=J,ncol=K)
 for(j in 1:J){
  for(k in 1:K){
   for(n in 1:J){
    a <- 1
    for(t in 1:T){
     if(exp(psitk[t,k]*(logitB[n,t]-logitB[j,t]))==Inf) {
      a <- a*1e6
    }else{ a <- a*exp(psitk[t,k]*(logitB[n,t]-logitB[j,t])) }  }
    B[j,k] <- B[j,k] + a  }
   B[j,k] <- 1/B[j,k] }  }

 B <- t(t(B)/colSums(B))
 # B is a JxK matrix
 pimais <- rowSums(obj)/sum(obj)

 P <- obj/rowSums(obj)   

 pk <- pimais %*% A # budget proportions 

 colnames(A) <- colnames(B) <- colnames(psitk) <- colnames(pk) <- paste('LB',1:K,sep='') 

 rownames(A) <- rownames(P)
 rownames(B) <- colnames(P)

 pij <- A %*% t(B) # expected budget

 residual <- P - pij

 val_func <- xab$value

 iter_ide <- round(as.numeric(xab$counts[2])/xab$outer.iterations) + xab$outer.iterations

 rescB <- rescaleB(obj,
                   A,
                   B)

 colnames(rescB) <- colnames(B)
 rownames(rescB) <- rownames(B)

 results <- list(P,
                 pij,
                 residual,
                 A, 
                 B,
                 rescB,
                 pk,
                 val_func,
                 iter_ide,
                 psitk)

 names(results) <- c('P',
                     'pij',
                     'residual',
                     'A',
                     'B',
                     'rescB',
                     'pk',
                     'val_func',
                     'iter_ide',
                     'psitk')

 #  names(results) <- c('Composition data matrix',
 #                      'Expected budget',
 #                      'Residual matrix',
 #                      'Mixing parameters',
 #                      'Latent budgets',
 #                      'Budget proportions',
 #                      'Value of the ls function',
 #                      'Number of iteractions',
 #                      'Logit parameters psitk')
 # 
 class(results)  <- c("lba.mle.logit",
                      "lba.mle")
 invisible(results)  
} 

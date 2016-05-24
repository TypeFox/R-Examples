lba.ls.fe <- function(obj          ,  
                      A            ,
                      B            ,
                      K            ,            
                      cA           ,
                      cB           ,
                      row.weights  ,
                      col.weights  ,
                      itmax.ide    ,
                      trace.lba    ,
                      ...)
{

 I  <- nrow(obj)       # row numbers of data matrix
 J  <- ncol(obj)       # column numbers of data matrix
 P  <- obj/rowSums(obj)  # observed components I x J

 #BUILDING aki
 if(!is.null(cA)& is.null(A)){  #if caki is given, otherwise only cbjk is given.
  #Find the indices of the matrix caki which are NA. Those indices will
  #be the ones, in which the matrix aki will be generated freely, the others
  #will be fixed or equal among them.
  A <- constrainAB(cA) }
 if(is.null(cA) & is.null(A)){
  #creating random generated values for alpha(k|i)   
  A <- rdirich(I, runif(K)) }

 #BUILDING bjk
 if(!is.null(cB) & is.null(B)){
  B <- t(constrainAB(t(cB))) }
 if(is.null(cB) & is.null(B)){
  #creating random generated values for beta(j|k) 
  B <- t(rdirich(K, runif(J))) }

 if(!is.null(cA)){
  minca <- min(cA, 
               na.rm=TRUE)
  maxca <- max(cA, 
               na.rm=TRUE) }
 if(!is.null(cB)) {
  mincb <- min(cB, 
               na.rm=TRUE)
  maxcb <- max(cB, 
               na.rm=TRUE) }

 #============================================================================    

 #============================================================================ 
 wls <- function(x, 
                 cA, 
                 cB, 
                 P, 
                 K, 
                 I, 
                 J, 
                 row.weights, 
                 col.weights ){
  y <- length(x)- (I*K+J*K)
  A <- matrix(x[(y+J*K +1):(y+J*K+I*K)], 
              ncol = K)
  B <- matrix(x[(y+1):(y+J*K)], 
              ncol = K)

  # Generating the identity matrix if row.weights is not informed
  if(is.null(row.weights)){
   V  <- diag(I)
  } else {
   V <- row.weights * diag(I)}

  # Generating the identity matrix if col.weights is not informed 
  if(is.null(col.weights)){
   W  <- diag(J)
  } else {
   wi <- col.weights
   W <- col.weights * diag(J) }

  ab <- A%*%t(B)
  wls <-  sum((V%*%(P - ab)%*%W)^2)
 }

 #============================================================================
 #         heq function
 #============================================================================
 heq <- function(x, 
                 cA, 
                 cB, 
                 P, 
                 K, 
                 I,
                 J,
                 row.weights,
                 col.weights) {

  # construction of matrices A(aki) and B(bjk) from input vector x
  y <- length(x)- (I * K + J * K )
  A <- matrix(x[(y+J*K +1):(y+J*K+I*K)],
              ncol = K)
  B <- matrix(x[(y+1):(y+J*K)],
              ncol = K)

  #constraints on bjk (B)
  if(!is.null(cB)) {
   mincb <- min(cB,
                na.rm=TRUE)
   maxcb <- max(cB,
                na.rm=TRUE) 
   if(maxcb > 1){                                       
    bl <- list()
    for(i in 2:max(cB, 
                   na.rm=TRUE)){
     bl[[i-1]] <- which(cB==i, 
                        arr.ind=TRUE)
     bl[[i-1]] <- bl[[i-1]][order(bl[[i-1]][,1],
                                  bl[[i-1]][,2]),]
    }  
    #this code forces the corresponding betasjk's to be equal
    h <- 0
    b <- 0
    e <- 0
    for(j in 1:(max(cB, 
                    na.rm=TRUE)-1)) {

     b[j] <- nrow(bl[[j]])

     for(i in 1:(b[j])){  e <- e+1

     h[e]<- B[bl[[j]][i,1],
              bl[[j]][i,2]]- x[sum(b)] 

     } 
    }  
   } 

   #this code forces the fixed constraints to be preserved on bjk (B).                      
   if(mincb <= 1) {

    posFb   <- which(cB<1, arr.ind = T)

    if(maxcb > 1){

     h[(length(h)+1):(length(h)+ nrow(posFb))] <- (B[posFb] - cB[posFb]) 

    }else{

     h <- 0  
     h[1:nrow(posFb)] <- (B[posFb] - cB[posFb]) 

    } 
   }  
  } 

  #constraints on aki (A)
  if(!is.null(cA)){
   minca <- min(cA, 
                na.rm=TRUE)
   maxca <- max(cA, 
                na.rm=TRUE) 
   if(maxca > 1){

    al <- list()

    for(i in 2:max(cA, 
                   na.rm=TRUE)){

     al[[i-1]] <- which(cA==i, 
                        arr.ind=TRUE)
     al[[i-1]] <- al[[i-1]][order(al[[i-1]][,1],
                                  al[[i-1]][,2]),]} 
    if(!is.null(cB)){ 

     if(maxcb > 1)  v <- length(h)                      

     if(maxcb <= 1)  h <- b <- v <- 0 

    }else{ h <- b <- v <- 0 }
    #this code forces the corresponding alphasik's to be equal
    a <- 0
    e <- 0

    for(j in 1:(max(cA, 
                    na.rm=TRUE)-1)) {

     a[j] <- nrow(al[[j]])

     for(i in 1:(a[j])){ e <- e+1

     h[v+e]<- A[al[[j]][i,1],al[[j]][i,2]]- x[sum(b)+sum(a)]

     }
    }  

    #this code forces the fixed constraints to be preserved on aki (A)
    if(minca <= 1){

     posFa   <- which(cA<1, arr.ind = T)

     h[(length(h)+1):(length(h)+ nrow(posFa))] <- (A[posFa] - cA[posFa])

    }
   }else{

    posFa   <- which(cA<1, arr.ind = T)

    if(!is.null(cB)) {

     h[(length(h)+1):(length(h)+ nrow(posFa))] <- (A[posFa] - cA[posFa])

    }else{

     h <- 0
     h[1:nrow(posFa)] <-  (A[posFa] - cA[posFa])

    } 
   } 
  } 

  #no constraints       
  if(all(c(is.null(cA),is.null(cB)))){ 

   h <- c(rowSums(A),colSums(B)) - rep(1,(I+K))
   #constraints    
  }else{

   h[(length(h)+1):(length(h)+(I+K))] <- c(rowSums(A),colSums(B)) - rep(1,(I+K))} 
  h
 }

 #=========================================================================
 #                hin  function
 #=========================================================================
 hin <- function(x, 
                 cA, 
                 cB, 
                 P, 
                 K, 
                 I, 
                 J, 
                 row.weights, 
                 col.weights) {
  y <- length(x)-(I*K+J*K)
  h <- x[(y+1):(y+I*K+J*K)] + 1e-9
  h
 }

 #===========================================================================
 #===========================================================================
 if(!is.null(cA)){
  if(maxca > 1){ #there are equality parameters in caki 

   #list containing in each element the positions of the equality parameters
   #of matrix A that are equal among them.
   al <- list()
   for(i in 2:max(cA, 
                  na.rm=TRUE)){
    al[[i-1]] <- which(cA==i, 
                       arr.ind=TRUE)
    al[[i-1]] <- al[[i-1]][order(al[[i-1]][,1],
                                 al[[i-1]][,2]),] }

   ma <- sum(sapply(al, function(x) nrow(x))) 

  } else { ma <- 0 } 
 } else { ma <- 0}

 if(!is.null(cB)){

  if(maxcb > 1) { #there are equality parameters in cbjk
   #list containing in each element the positions of the equality parameters
   #of matrix B that are equal among them.
   bl <- list()

   for(i in 2:max(cB, na.rm=TRUE)){

    bl[[i-1]] <- which(cB==i, 
                       arr.ind=TRUE)
    bl[[i-1]] <- bl[[i-1]][order(bl[[i-1]][,1],
                                 bl[[i-1]][,2]),]
   }

   mb <- sum(sapply(bl, 
                    function(x) nrow(x)))

  } else { mb <- 0 } 
 } else { mb <- 0 } 

 m <- ma + mb
 a <- rep(1e-6,m)
 x0 <- c(a,as.vector(B),
         as.vector(A))

 # finding the ls estimates
 itmax.ala <- round(.1*itmax.ide)
 itmax.opt <- round(.9*itmax.ide)
  
 xab <- constrOptim.nl(par     = x0,
                       fn      = wls,
                       cA      = cA,
                       cB      = cB,
                       P       = P,
                       K       = K,
                       I       = I,
                       J       = J,
                       row.weights  = row.weights,
                       col.weights  = col.weights,          
                       heq     = heq,
                       hin     = hin,
                       control.outer=list(trace=trace.lba,
                                          itmax=itmax.ala),
                       control.optim=list(maxit=itmax.opt))

 y <- length(xab$par)- (I * K + J * K )
 A <- matrix(xab$par[(y+J*K +1):(y+J*K+I*K)], 
             ncol = K)
 B <- matrix(xab$par[(y+1):(y+J*K)], 
             ncol = K)

 rownames(A) <- rownames(P)
 rownames(B) <- colnames(P)

 colnames(A) <- colnames(B) <- paste('LB',
                                     1:K,
                                     sep='')

 pimais <- rowSums(obj)/sum(obj)

 pk <- pimais %*% A # budget proportions

 colnames(pk) <- paste('LB',
                       1:K,
                       sep='')

 pij <- A %*% t(B) # expected budget

rownames(pij) <- rownames(P)
colnames(pij) <- colnames(P)

 residual <- P - pij

 val_func <- xab$value

 iter_ide <- round(as.numeric(xab$counts[2]/xab$outer.iterations)) + as.numeric(xab$outer.iterations)

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
                 iter_ide)

 names(results) <- c('P',
                     'pij',
                     'residual',
                     'A',
                     'B',
                     'rescB',
                     'pk',
                     'val_func',
                     'iter_ide') 

 class(results) <- c("lba.ls.fe",
                     "lba.ls")

 invisible(results)

}

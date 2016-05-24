lba.mle.fe <- function(obj      ,
                       A        ,
                       B        ,
                       K        , 
                       cA       ,
                       cB       ,
                       tolG     ,
                       tolA     ,
                       tolB     ,
                       itmax.ide,
                       trace.lba,
                       toltype  ,  
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
 #The default toltype is "all".
 I <- nrow(obj)
 J <- ncol(obj)

 P <- obj/rowSums(obj)
 #-----------------------------------------------------------------------------

 #BUILDING aki (A)
 if(!is.null(cA) & is.null(A)){  
  #Find the indices of the matrix caki which are NA. Those indices will
  #be the ones, in which the matrix aki (A) will be generated freely, the others
  #will be fixed or equal among them.
  A <- constrainAB(cA) 
 }

 if(is.null(cA) & is.null(A)){
  #creating random generated values for alpha(i|k) using Dirichlet function
  A <- rdirich(I, runif(K))  
 }

 #BUILDING bjk (B)

 if(all(!is.null(cB), is.null(B))){
  #Find the indices of the matrix cbjk which are NA. Those indices will
  #be the ones, in which the matrix bjk (B) will be generated freely, the others
  #will be fixed or equal among them.
  B <- t(constrainAB(t(cB)))  
 }
 
 if(all(c(is.null(cB), is.null(B)))){
  #creating random generated values for beta(j|k) using Dirichlet function
  B <- t(rdirich(K, runif(J))) 
 }

 if(!is.null(cA)){
  maxca <- max(cA, na.rm=TRUE)
  minca <- min(cA, na.rm=TRUE)
  ca <- cA
  ca[ca>1] <- NA  
 }    #taking out the equality parameters 

 if(!is.null(cB)){
  maxcb <- max(cB, na.rm=TRUE)
  mincb <- min(cB, na.rm=TRUE)
  cb <- cB
  cb[cb>1] <- NA 
 }   #taking out the equality parameters 
 
 #-----------------------------------------------------------------------------
 #starting MLE algorithm
 iter_ide <- 0
 
 repeat
 {
  iter_ide <- iter_ide + 1
  #-------------------------------------------------------------------------------
  #previous G2 to be used in |G2 - G2a|<tol
  qij <- A%*%t(B)
  kij <- obj/rowSums(obj)
  G2a <- 2*sum(obj*log(kij/qij))
  akia <- A
  bjka <- B
  #-------------------------------------------------------------------------------
  ini <- list()
  #----------------------------------------------------------------------------
  #for(k in 1:K) ini[[k]] <- c(a[,k],b[,k])
  #-------------------------------------------------------------------------------
  for(k in 1:K) ini[[k]] <- c(A[,k],B[,k])
  #------------------------------------------------------------------------------
  #initial values list with vectors a1(I),b1(J);... ak(I),bk(J);
  # ... aK(I), bK(J) where a's are alfa values and b's are beta values. The lenghts
  # of the a's are I and of the b's are J. The elements of the list are c(a1,b1),
  #... c(ak,bk),... c(aK,bK).
  #
  pijk <-   lapply(ini, function(xx) outer(xx[1:I],xx[(I+1):(I+J)])*(rowSums(obj)/sum(obj)))
  #The elements of the list pijk are matrices (I x J) for fixed k
  #and the number of elements is K. They are alpha(k|i)beta(j|k)*pi+ =
  #pijk|i*pi+
  #
  mpijk <- sapply(pijk, function(x) x, simplify=TRUE)
  #Transform pijk into a matrix I*JxK. The columns are the colmuns of each element
  #(matrix) one after the other, of the first then the second and so on.
  #
  sij <- rowSums(mpijk)
  #The sum over k of alpha(k|i)beta(j|k), vector of lenght IxJ
  #
  lnijk <- lapply(pijk, function(xx) obj*xx )
  #The product of N(ij)*pijk. It is a list of K matrices IxJ each
  #
  nijk <- lapply(lnijk, function(x) x/sij )
  #The list containing the values of n(hatzero)(ijk). Each element is a matrix IxJ and
  #the number of elements of the list is K.
  #
  sjnijk <- sapply(nijk, function(x) rowSums(x), simplify=TRUE)
  #sum over j for i and k fixed. sjnijk is a IxK matrix. ni+k
  #
  sinijk <- sapply(nijk, function(x) colSums(x), simplify=TRUE)
  #sum over i for j and k fixed. sinijk is a JxK matrix  n+jk

  #------------------------------------------------------------------------------
  #CONSTRAINTS ON caki (A) 
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------    
  if(!is.null(cA)){
   #------------------------------------------------------------------------------
   # ONLY FIXED CONSTRAINTS 
   if(maxca<=1){

    A <- constmleFA(cA,
                    A,
                    sjnijk)

   }  
   # next values alphahat1 to be used in the next step of the iteration

   #------------------------------------------------------------------------------
   # ONLY EQUALITY CONSTRAINTS 
   if(minca > 1){

    A <- constmleEA(cA,
                    A,
                    sjnijk)

   } 
   # next values alphahat1 to be used in the next step of the iteration
   #------------------------------------------------------------------------------

   #------------------------------------------------------------------------------
   # FIXED AND EQUALITY CONSTRAINTS
   if(all(c(minca<=1, maxca>1))) {
    #list containing in each element the positions of the equality parameters
    #of matrix A that are equal among them.
    al <- list()

    for(i in 2:max(cA, na.rm=TRUE)){
     al[[i-1]] <- which(cA==i, arr.ind=TRUE)
     al[[i-1]] <- al[[i-1]][order(al[[i-1]][,1],al[[i-1]][,2]),] 
    }

    if(is.null(cB)){
     if(any(all(sapply(al,function(x) all(x[,2] == x[1,2]))),
            all(ca[!is.na(ca)] == 0))){

      A <- constmleFEA(cA,
                       A,
                       sjnijk)

      # next values alphahat1 to be used in the next step of the iteration 
     } else {
      res <- constmleFEalabama(obj, 
                               A, 
                               B, 
                               cA, 
                               cB,
                               itmax.ide,
                               trace.lba)
      return(res)
     }  
    }else{
     if(mincb>1){
      if(any(all(sapply(al,function(x) all(x[,2] == x[1,2]))),
             all(ca[!is.na(ca)] == 0))){

       A <- constmleFEA(cA,
                        A,
                        sjnijk)


       # next values alphahat1 to be used in the next step of the iteration 
      } else {

       res <- constmleFEalabama(obj, 
                                A, 
                                B, 
                                cA,
                                cB,
                                itmax.ide,
                                trace.lba)
       return(res)
      }  
     }

     if(maxcb<=1){
      if(any(all(sapply(al,function(x) all(x[,2] == x[1,2]))),
             all(ca[!is.na(ca)] == 0))){

       A <- constmleFEA(cA,
                        A,
                        sjnijk)

       # next values alphahat1 to be used in the next step of the iteration 
      } else {
       res <- constmleFEalabama(obj, 
                                A, 
                                B, 
                                cA, 
                                cB,
                                itmax.ide,
                                trace.lba)
       return(res)
      } 
     }    

     if(all(c(mincb<=1, maxcb>1))){ 
      bl <- list()
      for(i in 2:max(cB, na.rm=TRUE)){
       bl[[i-1]] <- which(cB==i, arr.ind=TRUE)
       bl[[i-1]] <- bl[[i-1]][order(bl[[i-1]][,1],bl[[i-1]][,2]),]
      } 

      if(any(all(sapply(bl,function(x) all(x[,2] == x[1,2]))),
             all(cb[!is.na(cb)] == 0))){   

       if(any(all(sapply(al,function(x) all(x[,2] == x[1,2]))),
              all(ca[!is.na(ca)] == 0))){

        A <- constmleFEA(cA,
                         A,
                         sjnijk)

        # next values alphahat1 to be used in the next step of the iteration 
       } else {
        res <- constmleFEalabama(obj, 
                                 A, 
                                 B, 
                                 cA, 
                                 cB,
                                 itmax.ide,
                                 trace.lba)
        return(res)
       }  
      }  
     } 
    } 
   }
  }
  #------------------------------------------------------------------------------     
  #------------------------------------------------------------------------------
  #CONSTRAINTS ON cbjk (B) 
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------ 
  if(!is.null(cB)){

   #------------------------------------------------------------------------------
   # ONLY FIXED CONSTRAINTS 
   if(maxcb<=1) {

    B <- constmleFB(cB,
                    B,
                    sinijk)

   }  

   # next values betahat1 to be used in the next step of the iteration
   #------------------------------------------------------------------------------
   # ONLY EQUALITY CONSTRAINTS 
   if(mincb > 1){

    B <- constmleEB(cB,
                    B,
                    sinijk)

   }  
   # next values betahat1 to be used in the next step of the iteration

   #------------------------------------------------------------------------------
   # FIXED AND EQUALITY CONSTRAINTS 
   if(all(c(mincb<=1, maxcb>1))){ 
    #list containing in each element the positions of the equality parameters
    #of matrix B that are equal among them.
    bl <- list()
    for(i in 2:max(cB, na.rm=TRUE)){
    
     bl[[i-1]] <- which(cB==i, arr.ind=TRUE)
     bl[[i-1]] <- bl[[i-1]][order(bl[[i-1]][,1],bl[[i-1]][,2]),]
    
    }
   
    if(is.null(cA)){

     if(any(all(sapply(bl,function(x) all(x[,2] == x[1,2]))),
            all(cb[!is.na(cb)] == 0))){

      B <- constmleFEB(cB,
                       B,
                       sinijk)

      # next values betahat1 to be used in the next step of the iteration
     } else {
      res <- constmleFEalabama(obj, 
                               A, 
                               B, 
                               cA, 
                               cB,
                               itmax.ide,
                               trace.lba)
      return(res)

     }     
    } else {

     if(mincb>1){
      if(any(all(sapply(bl,function(x) all(x[,2] == x[1,2]))),
             all(cb[!is.na(cb)] == 0))){

       B <- constmleFEB(cB,
                        B,
                        sinijk)

       # next values betahat1 to be used in the next step of the iteration
      } else {
       res <- constmleFEalabama(obj, 
                                A, 
                                B, 
                                cA, 
                                cB,
                                itmax.ide,
                                trace.lba)
       return(res)
      } 
     }

     if(maxca<=1){
      if(any(all(sapply(bl,function(x) all(x[,2] == x[1,2]))),
             all(cb[!is.na(cb)] == 0))){

       B <- constmleFEB(cB,
                        B,
                        sinijk)

       # next values betahat1 to be used in the next step of the iteration
      } else {
       res <- constmleFEalabama(obj, 
                                A, 
                                B, 
                                cA, 
                                cB,
                                itmax.ide,
                                trace.lba)
       return(res)
      } 
     }

     if(all(c(minca<=1, maxca>1))) { 
      
      al <- list()
      
      for(i in 2:max(cA, na.rm=TRUE)){
      
       al[[i-1]] <- which(cA==i, arr.ind=TRUE)
       al[[i-1]] <- al[[i-1]][order(al[[i-1]][,1],al[[i-1]][,2]),] 
      
      } 

      if(any(all(sapply(al,function(x) all(x[,2] == x[1,2]))),
             all(ca[!is.na(ca)] == 0))){ 

       if(any(all(sapply(bl,function(x) all(x[,2] == x[1,2]))),
              all(cb[!is.na(cb)] == 0))){

        B <- constmleFEB(cB,
                         B,
                         sinijk)

        # next values betahat1 to be used in the next step of the iteration
       } else {
        res <- constmleFEalabama(obj, 
                                 A, 
                                 B, 
                                 cA, 
                                 cB,
                                 itmax.ide,
                                 trace.lba)
        return(res)
       } 
      } 
     } 
    } 
   } 
  # FIXED AND EQUALITY CONSTRAINTS ALABAMA
  #=============================================================================
  if(all(!is.null(cA),!is.null(cB))){ 
   if(all(minca<=1, maxca>1)){
    
    al <- list()
    
    for(i in 2:max(cA, na.rm=TRUE)){
    
     al[[i-1]] <- which(cA==i, arr.ind=TRUE)
     al[[i-1]] <- al[[i-1]][order(al[[i-1]][,1],al[[i-1]][,2]),] 
    
    }

    if(any(sapply(al,function(x) any(x[,2] != x[1,2])))){                        

     if(all(mincb<=1, maxcb>1)){
     
      bl <- list()
      
      for(i in 2:max(cB, na.rm=TRUE)){
      
       bl[[i-1]] <- which(cB==i, arr.ind=TRUE)
       bl[[i-1]] <- bl[[i-1]][order(bl[[i-1]][,1],bl[[i-1]][,2]),]
      
      }  

      if(any(sapply(bl,function(x) any(x[,2] != x[1,2])))){                          

       res <- constmleFEalabama(obj, 
                                A, 
                                B, 
                                cA, 
                                cB,
                                itmax.ide,
                                trace.lba)
       return(res)
      } 
     } 
    } 
   } 
  } 
  }
  #-----------------------------------------------------------------------------

  pij <- A%*%t(B)
  #Values of pihat(j|i) = the sum over k of alphahat(k|i)*betahat(j|k)
  G2 <- 2*sum(obj*log(obj/(pij*rowSums(obj))))
  chi2 <- sum(((obj - pij*rowSums(obj))^2)/obj)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  at <- max(abs(A - akia))
  bt <- max(abs(B - bjka))
  Gt <- abs(G2 - G2a)
  #------------------------------------------------------------------------------
  if(iter_ide > itmax.ide){warning("maximum number of iteractions exceeded" )
  break }

  if (toltype == "all" & Gt < tolG & at < tolA & bt < tolB) {
   break  #default
  } 

  if (toltype == "G2" & Gt < tolG) {
    break
   } 

  if (toltype == "ab" & at < tolA & bt < tolB) {
     break
    }
   
 #------------------------------------------------------------------------------
 #closing the repeat
}

pip <- rowSums(obj)/sum(obj) #this is pi+
pjki <- rep(0,I*J*K) #this will become pijk/pi+ see van der Ark page 80
m <- 0
# this makes i=1,j=1,k=1; i=2,j=1,k=1; i=3,j=1,k=1...; i=2,j=1,k=1 and so on.
for(k in 1:K) for(j in 1:J) for(i in 1:I) { m <-  m +1
pjki[m] <- A[i,k]*B[j,k] }
pjki[pjki==0] <- 1e-7
mi <- matrix(pjki, I)
pijk <- mi*pip #this is pijk see van der Ark page 80
nijk <- as.vector(pijk*sum(obj))
val_func <- -sum(nijk * log(pjki)) # -loglikelihood function

pij <- A %*% t(B)

residual <- P - pij

pk <- pip %*% A  # budget proportions

rownames(A) <- rownames(P)
rownames(B) <- colnames(P)

colnames(A) <- colnames(B) <- colnames(pk) <- paste('LB',1:K,sep='')

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

#  names(results) <- c('Composition data matrix',
#                      'Expected budget',
#                      'Residual matrix',
#                      'Mixing parameters',
#                      'Latent budgets',
#                      'Budget proportions',
#                      'Value of the -loglik function',
#                      'Number of iteractions')
# 
class(results) <- "lba.mle.fe"

invisible(results)

}

constmleFA <- function(cA, 
                       A, 
                       sjnijk)
{
 awf <- which(is.na(cA), arr.ind=TRUE)
 awf <- awf[ order(awf[,1],awf[,2]), ] #indices of caki = NA

 ac <- sjnijk[awf]
 cjknijk <- numeric()
 for(i in 1:nrow(cA)) cjknijk[i] <- sum(ac[awf[,1]==i])
 #sum over j and k and i fixed. cjknijk is a vector of lenght I, only free parameters 

 A[awf] <- (sjnijk/cjknijk)[awf]

 for(i in 1:nrow(cA))  if(sum(A[i,][!is.na(cA[i,]) < 1])){
  sba <- sum(A[i,][is.na(cA[i,])])
  sba.1 <- 1 - sum(A[i,][!is.na(cA[i,])])
  A[i,][is.na(cA[i,])] <- A[i,][is.na(cA[i,])]/sba*sba.1 

 }
 #the values of aki[awf] were adjusted to satisfy the constraints sum = 1
 #now the aki[awf] are next values alphahat1 to be used in the next step of the 
 #iteration. aki[awf] corresponds only to the free parameters.
 invisible(A)
}  

constmleFB <- function(cB, 
                       B, 
                       sinijk)
{

 bwf <- which(is.na(t(cB)), arr.ind=TRUE)
 bwf <- bwf[ order(bwf[,1],bwf[,2]), ] #indices of cbjk = NA

 bc <- t(sinijk)[bwf]
 cijnijk <- numeric()
 for(i in 1:ncol(cB)) cijnijk[i] <- sum(bc[bwf[,1]==i])
 #sum over i and j and k fixed. cijnijk is a vector of lenght K, only free parameters

 bsc <- t(sinijk)/cijnijk
 tbjk <- t(B)
 tbjk[bwf] <- bsc[bwf]
 B <- t(tbjk)
 for(i in 1:ncol(cB))  if(sum(B[,i][!is.na(cB[,i]) < 1])){
  sab <- sum(B[,i][is.na(cB[,i])])
  sab.1 <- 1 - sum(B[,i][!is.na(cB[,i])])
  B[,i][is.na(cB[,i])] <- B[,i][is.na(cB[,i])]/sab*sab.1}

 #the values of bjk[bwf] were adjusted to satisfy the constraints sum = 1
 #now the bjk[bwf] are next values beta hat1 to be used in the next step of the
 #iteration. bjk[bwf] corresponds only to the free parameters.
 invisible(B)
} 

constmleEA <- function(cA, 
                       A, 
                       sjnijk)
{

 awe <- which(is.na(cA), arr.ind=TRUE)
 awe <- awe[ order(awe[,1],awe[,2]), ] #positions of NA's

 sjknijk <- apply(sjnijk, 1, function(x) sum(x))
 # sum over j and k and i fixed. sjknijk is a vector of lenght I. ni++
 al <- list()
 a1 <- cA
 for(i in 2:max(a1, na.rm=TRUE)){
  al[[i-1]] <- which(a1==i, arr.ind=TRUE)
  al[[i-1]] <- al[[i-1]][
                         order(al[[i-1]][,1],al[[i-1]][,2]),]
 }
 #list containing in each element the positions of the equality parameters
 #that are equal among them.
 #

 sapply(al, function(x) A[x]<<- sum( sjnijk[x])/sum(sjknijk[x[,1]]))
 #the update of the equality parameters to be used in the next step
 #
 A[awe] <- (sjnijk/sjknijk)[awe]

 if(all(is.na(rowSums(cA)))){
  #only if all rows aki have at least one NA

  for(i in 1:nrow(cA)) A[i,] <- A[i,]/rowSums(A)[i]
  for(i in 1:nrow(cA))  if(sum(A[i,][!is.na(cA[i,]) < 1])){
   sba <- sum(A[i,][is.na(cA[i,])])
   sba.1 <- 1 - sum(A[i,][!is.na(cA[i,])])
   A[i,][is.na(cA[i,])] <-
    A[i,][is.na(cA[i,])]/sba*sba.1 

  }
  #adjusted random values to satisfy the constraints sum = 1
  # next values alphahat1 to be used in the next step of the iteration
  #of the free parameters.
 }else{
  # if aki has at least one row without NA´s

  ra <- which(!is.na(rowSums(cA)))
  a <- A

  for(i in ra) A[i,] <- A[i,]/rowSums(A)[i]
  for(i in ra) for(j in 1:ncol(cA)) A[a==a[i,j]] <- A[i,j]

  sapply(al, function(x) A[x] <<- min(A[x]))
  for(i in 1:nrow(cA))  if(sum(A[i,][!is.na(cA[i,]) < 1])){

   sba <- sum(A[i,][is.na(cA[i,])])
   sba.1 <- 1 - sum(A[i,][!is.na(cA[i,])])
   A[i,][is.na(cA[i,])] <-
    A[i,][is.na(cA[i,])]/sba*sba.1 

  }

  #adjusted random values to satisfy the constraints sum = 1
  # next values alphahat1 to be used in the next step of the iteration
  #of the free parameters.
  if(any(A < 0)) {
   A <- constrainAB(cA)
  }
 }
 invisible(A)
} 

constmleEB <- function(cB, 
                       B, 
                       sinijk)
{

 bwe <- which(is.na(cB), arr.ind=TRUE)
 bwe <- bwe[ order(bwe[,1],bwe[,2]), ] #positions of NA's

 sijnijk <- apply(sinijk, 2, function(x) sum(x))
 # sum over i and j and k fixed. sjknijk is a vector of lenght K.
 bl <- list()
 b1 <- cB

 for(i in 2:max(b1, na.rm=TRUE)){
  bl[[i-1]] <- which(b1==i, arr.ind=TRUE)
  bl[[i-1]] <- bl[[i-1]][
                         order(bl[[i-1]][,1],bl[[i-1]][,2]),]
 }
 #list containing in each element the positions of the equality parameters
 #that are equal among them.

 sapply(bl, function(x) B[x]<<- sum( sinijk[x])/sum(sijnijk[x[,2]]))
 #the update of the equality parameters to be used in the next step

 B[bwe] <- t(t(sinijk)/sijnijk)[bwe]

 if(all(is.na(colSums(cB)))){
  #only if all colmuns of cbjk have at least one NA

  for(i in 1:ncol(cB)) B[,i] <- B[,i]/colSums(B)[i]
  sapply(bl, function(x) B[x] <<- min(B[x]))

  for(i in 1:ncol(cB))  if(sum(B[,i][!is.na(cB[,i]) < 1])){

   sba <- sum(B[,i][is.na(cB[,i])])
   sba.1 <- 1 - sum(B[,i][!is.na(cB[,i])])
   B[,i][is.na(cB[,i])] <- B[,i][is.na(cB[,i])]/sba*sba.1 

  }
  # next values betahat1 to be used in the next step of the iteration
  #of the free parameters.

 }else{
  # if bjk has at least one column without NA´s

  rb <- which(!is.na(colSums(cB)))

  b <- B

  for(i in rb) B[,i] <- B[,i]/sum(B[,i])
  for(i in rb) for(j in 1:nrow(cB)) B[b==b[j,i]] <- B[j,i]

  for(i in 1:ncol(cB))  if(sum(B[,i][!is.na(cB[,i]) < 1])){
   sba <- sum(B[,i][is.na(cB[,i])])
   sba.1 <- 1 - sum(B[,i][!is.na(cB[,i])])
   B[,i][is.na(cB[,i])] <- B[,i][is.na(cB[,i])]/sba*sba.1 }
  #adjusted random values to satisfy the constraints sum = 1
  # next values alphahat1 to be used in the next step of the iteration
  #of the free parameters.
  if(any(B < 0)){

   B <- t(constrainAB(t(cB)))

  }
 }
 invisible(B)
} 

constmleFEA <- function(cA, 
                        A, 
                        sjnijk)
{
 K <- ncol(cA)
 ca <- cA
 ca[ca>1] <- NA      #taking out the equality parameters
 awf <- which(is.na(ca), arr.ind=TRUE)
 awf <- awf[ order(awf[,1],awf[,2]), ] #indices of ca = NA
 awe <- which(is.na(cA), arr.ind=TRUE) # back the equality parameters
 awe <- awe[ order(awe[,1],awe[,2]), ] #positions of NA's
 posFa   <- which(cA<1, arr.ind = T)

 ac <- sjnijk[awf]
 cjknijk <- numeric()
 for(i in 1:nrow(cA)) cjknijk[i] <- sum(ac[awf[,1]==i])
 #sum over j and k and i fixed. cjknijk is a vector of lenght I, only free parameters

 A[awf] <- (sjnijk/cjknijk)[awf]
 #only equality and free parameters
 #-----------------------------------------------------------------------------
 al <- list()
 a1 <- cA
 for(i in 2:max(a1, na.rm=TRUE)){
  al[[i-1]] <- which(a1==i, arr.ind=TRUE)
  al[[i-1]] <- al[[i-1]][order(al[[i-1]][,1],al[[i-1]][,2]),]
 }
 #list containing in each element the positions of the equality parameters
 #that are equal among them.

 lapply(al, function(x) A[x]<<- sum( sjnijk[x])/sum(cjknijk[x[,1]]))
 #the update of the equality parameters to be used in the next step

 A[awe] <- (sjnijk/cjknijk)[awe]
 #    next values alphahat1 to be used in the next step of the iteration
 #   of the free parameters.

 if(all(is.na(rowSums(cA)))){
  #only if all rows caki have at least one NA

  for(i in 1:nrow(cA)){
   A[i,][is.na(ca[i,])] <- A[i,][is.na(ca[i,])]* (1-sum(A[i,][!is.na(ca[i,])]))/(sum(A[i,][is.na(ca[i,])])) 
  }

  sapply(al, function(x) A[x] <<- min(A[x]))

  for(i in 1:nrow(cA))  if(sum(A[i,][!is.na(cA[i,]) < 1])){
   sba <- sum(A[i,][is.na(cA[i,])])
   sba.1 <- 1 - sum(A[i,][!is.na(cA[i,])])
   A[i,][is.na(cA[i,])] <-
    A[i,][is.na(cA[i,])]/sba*sba.1 }
  #adjusted random values to satisfy the constraints sum = 1
  # next values alphahat1 to be used in the next step of the iteration
  #of the free parameters.

 }else{
  # if aki has at least one row without NA´s

  ra <- which(!is.na(rowSums(cA)))

  a <- A
  for(i in ra) A[i,] <- A[i,]/sum(A[i,])
  for(i in ra) for(j in 1:K) A[a==a[i,j]] <- A[i,j]

  for(i in 1:nrow(cA)){ if(sum(A[i,][!is.na(cA[i,]) < 1])){
   sba <- sum(A[i,][is.na(cA[i,])])
   sba.1 <- 1 - sum(A[i,][!is.na(cA[i,])])
   A[i,][is.na(cA[i,])] <-
    A[i,][is.na(cA[i,])]/sba*sba.1 }
   #adjusted random values to satisfy the constraints sum = 1
   # next values alphahat1 to be used in the next step of the iteration
   #of the free parameters.

   if(any(A < 0)) {
    A <- constrainAB(cA)
   }
  }
 }
 invisible(A)
}

constmleFEB <- function(cB, 
                        B, 
                        sinijk)
{

 cb <- cB
 cb[cb>1] <- NA
 bwf <- which(is.na(t(cb)), arr.ind=TRUE)
 bwf <- bwf[ order(bwf[,1],bwf[,2]), ] #indices of cb = NA
 bwe <- which(is.na(cB), arr.ind=TRUE)
 bwe <- bwe[ order(bwe[,1],bwe[,2]), ] #positions of NA's
 posFb   <- which(cB<1, arr.ind = T)

 bc <- t(sinijk)[bwf]
 cijnijk <- numeric()
 for(i in 1:ncol(cB)) cijnijk[i] <- sum(bc[bwf[,1]==i])
 #sum over i and j and k fixed. cijnijk is a vector of lenght K, only free parameters

 bsc <- t(sinijk)/cijnijk
 tbjk <- t(B)
 tbjk[bwf] <- bsc[bwf]
 B <- t(tbjk)
 #only equality and free parameters

 bl <- list()
 b1 <- cB
 for(i in 2:max(b1, na.rm=TRUE)){
  bl[[i-1]] <- which(b1==i, arr.ind=TRUE)
  bl[[i-1]] <- bl[[i-1]][
                         order(bl[[i-1]][,1],bl[[i-1]][,2]),]
 }

 lapply(bl, function(x) B[x]<<- sum( sinijk[x])/sum(cijnijk[x[,2]]))
 #the update of the equality parameters to be used in the next step
 B[bwe] <- t(t(sinijk)/cijnijk)[bwe]
 # next values alphahat1 to be used in the next step of the iteration
 #of the free parameters.

 #    cjk <- cbjk
 #    cjk[posFb] <- NA
 if(all(is.na(colSums(cB)))){
  #only if all colmuns of cbjk have at least one NA

  for(i in 1:ncol(cB)){
   B[,i][is.na(cb[,i])] <- B[,i][is.na(cb[,i])]* (1-sum(B[,i][!is.na(cb[,i])]))/(sum( B[,i][is.na(cb[,i])])) 
  }

  sapply(bl, function(x) B[x] <<- min(B[x]))


  for(i in 1:ncol(cB))  if(sum(B[,i][!is.na(cB[,i]) < 1])){

   sba <- sum(B[,i][is.na(cB[,i])])
   sba.1 <- 1 - sum(B[,i][!is.na(cB[,i])])
   B[,i][is.na(cB[,i])] <- B[,i][is.na(cB[,i])]/sba*sba.1 

  }

  # next values betahat1 to be used in the next step of the iteration
  #of the free parameters.
 }else{
  # if bjk has at least one column without NA´s
  rb <- which(!is.na(colSums(cB)))
  b <- B
  for(i in rb) B[,i] <- B[,i]/sum(B[,i])
  for(j in rb) for(i in 1:nrow(cB)) B[b==b[i,j]] <- B[i,j]
  for(i in 1:ncol(cB)) if(sum(B[,i][!is.na(cB[,i]) < 1])){
   sba <- sum(B[,i][is.na(cB[,i])])
   sba.1 <- 1 - sum(B[,i][!is.na(cB[,i])])
   B[,i][is.na(cB[,i])] <- B[,i][is.na(cB[,i])]/sba*sba.1 }
  #adjusted random values to satisfy the constraints sum = 1
  # next values alphahat1 to be used in the next step of the iteration
  #of the free parameters.

  if(any(B < 0)){

   B <- t(constrainAB(t(cB)))

  }
 }
 invisible(B)
}

constmleFEalabama <- function(obj, 
                              A, 
                              B, 
                              cA, 
                              cB,
                              itmax.ide,
                              trace.lba)
{
 I <- nrow(A)
 J <- nrow(B)
 K <- ncol(A)
 P <- obj/rowSums(obj)
 if(!is.null(cA)) {
  minca <- min(cA, na.rm=TRUE)
  maxca <- max(cA, na.rm=TRUE) }
 if(!is.null(cB)) {  
  mincb <- min(cB, na.rm=TRUE)
  maxcb <- max(cB, na.rm=TRUE) }

 #=============================================================================
 # FIXED AND EQUALITY CONSTRAINTS ALABAMA
 #=============================================================================

 mle <- function(xx,obj,cA, cB, P, K, I, J){
  y <- length(xx)- (I*K+J*K)
  A <- matrix(xx[(y+J*K +1):(y+J*K+I*K)], ncol = K)
  B <- matrix(xx[(y+1):(y+J*K)], ncol = K)

  pip <- rowSums(obj)/sum(obj) #this is pi+
  pjki <- rep(0,I*J*K) #this will become pijk/pi+ see van der Ark page 80
  m <- 0
  # this makes i=1,j=1,k=1; i=2,j=1,k=1; i=3,j=1,k=1...; i=2,j=1,k=1 and so on.
  for(k in 1:K) for(j in 1:J) for(i in 1:I) { m <-  m +1
  pjki[m] <- A[i,k]*B[j,k] }
  pjki[pjki <= 0] <- 1e-7
  mi <- matrix(pjki, I)
  pijk <- mi*pip #this is pijk see van der Ark page 80
  nijk <- as.vector(pijk*sum(obj))
  mle <- -sum(nijk * log(pjki)) # -loglikelihood function
 }

 #============================================================================
 #         heq function
 #============================================================================
 heq <- function(xx,obj, cA, cB, P, K, I, J) {

  # construction of matrices A(A) and B(bjk) from input vector x
  y <- length(xx)- (I * K + J * K )
  A <- matrix(xx[(y+J*K +1):(y+J*K+I*K)], ncol = K)
  B <- matrix(xx[(y+1):(y+J*K)], ncol = K)

  #constraints on bjk (B)
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
    for(j in 1:(max(cB, na.rm=TRUE)-1)) {
     b[j] <- nrow(bl[[j]])
     for(i in 1:(b[j])){  e <- e+1
     h[e]<- B[bl[[j]][i,1],bl[[j]][i,2]]- xx[sum(b)] 
    } 
    }   

    #this code forces the fixed constraints to be preserved on bjk (B).                      
    if(mincb <= 1) {
     posFb   <- which(cB<1, arr.ind = T)
     h[(length(h)+1):(length(h)+ nrow(posFb))] <-(B[posFb] - cB[posFb])  
    }
   }else{
    posFb   <- which(cB<1, arr.ind = T)
    h <- 0  
    h[1:nrow(posFb)] <- (B[posFb] - cB[posFb]) 
  }  
  }  

  #constraints on A (A)
  if(!is.null(cA)){
   minca <- min(cA, na.rm=TRUE)
   maxca <- max(cA, na.rm=TRUE) 
   if(maxca > 1){
    al <- list()
    for(i in 2:max(cA, na.rm=TRUE)){
     al[[i-1]] <- which(cA==i, arr.ind=TRUE)
     al[[i-1]] <- al[[i-1]][
                            order(al[[i-1]][,1],al[[i-1]][,2]),]} 
    if(!is.null(cB)){ 
     if(maxcb > 1)  v <- length(h)                      
     if(maxcb <= 1)  h <- b <- v <- 0 
    }else{ h <- b <- v <- 0 }
    #this code forces the corresponding alphasik's to be equal
    a <- 0
    e <- 0
    for(j in 1:(max(cA, na.rm=TRUE)-1)) {
     a[j] <- nrow(al[[j]])
     for(i in 1:(a[j])){ e <- e+1
     h[v+e]<-
      A[al[[j]][i,1],al[[j]][i,2]]- xx[sum(b)+sum(a)] } }  

    #this code forces the fixed constraints to be preserved on A (A)
    if(minca <= 1){
     posFa   <- which(cA<1, arr.ind = T)
     h[(length(h)+1):(length(h)+ nrow(posFa))] <- (A[posFa] - cA[posFa]) 
    }
   }else{
    posFa   <- which(cA<1, arr.ind = T)
    if(!is.null(cB)) {
     h[(length(h)+1):(length(h)+ nrow(posFa))] <-(A[posFa] - cA[posFa])
    }else{
     h <- 0
     h[1:nrow(posFa)] <-  (A[posFa] - cA[posFa])
  } 
  } 
  }     

  h[(length(h)+1):(length(h)+(I+K))] <- c(rowSums(A),colSums(B)) - rep(1,(I+K))
  h
 }

 #=========================================================================
 #                hin  function
 #=========================================================================
 hin <- function(xx, obj, cA, cB, P, K, I, J) {
  y <- length(xx)-(I*K+J*K)
  h <- xx[(y+1):(y+I*K+J*K)] + 1e-9
  h
 }

 #===========================================================================
 #===========================================================================
 if(!is.null(cA)){
  if(maxca > 1){ #there are equality parameters in caki
   #list containing in each element the positions of the equality parameters
   #of matrix A that are equal among them.
   al <- list()
   for(i in 2:max(cA, na.rm=TRUE)){
    al[[i-1]] <- which(cA==i, arr.ind=TRUE)
    al[[i-1]] <- al[[i-1]][order(al[[i-1]][,1],al[[i-1]][,2]),] 
   }
   ma <- sum(sapply(al, function(x) nrow(x)))
  } 
 } else { ma <- 0 }                                     

 if(!is.null(cB)){
  if(maxcb > 1) { #there are equality parameters in cbjk
   #list containing in each element the positions of the equality parameters
   #of matrix B that are equal among them.
   bl <- list()
   for(i in 2:max(cB, na.rm=TRUE)){
    bl[[i-1]] <- which(cB==i, arr.ind=TRUE)
    bl[[i-1]] <- bl[[i-1]][order(bl[[i-1]][,1],bl[[i-1]][,2]),]
   }
   mb <- sum(sapply(bl, function(x) nrow(x))) 
  }
 } else { mb <- 0 }

 m <- ma + mb
 a <- rep(1e-6,m)
 x0 <- c(a,as.vector(B),as.vector(A))

 # finding the mle estimates
 itmax.ala <- round(0.1*itmax.ide)
 itmax.opt <- round(0.9*itmax.ide)
 # 
 xab <- constrOptim.nl(par     = x0,
                       fn      = mle,
                       cA      = cA,
                       cB      = cB,
                       obj     = obj,
                       P       = P,
                       K       = K,
                       I       = I,
                       J       = J,
                       heq     = heq,
                       hin     = hin,
                       control.outer=list(trace=trace.lba,
                                          itmax=itmax.ala),
                       control.optim=list(maxit=itmax.opt))
 # 
 y <- length(xab$par)- (I * K + J * K )
 A <- matrix(xab$par[(y+J*K +1):(y+J*K+I*K)], ncol = K)
 B <- matrix(xab$par[(y+1):(y+J*K)], ncol = K)

 iter_ide <- round(as.numeric(xab$counts[2])/xab$outer.iterations) + xab$outer.iterations

 pip <- rowSums(obj)/sum(obj) #this is pi+ 

 pij <- A %*% t(B)

 residual <- P - pij

 pk <- pip %*% A  # budget proportions

 rownames(A) <- rownames(P)
 rownames(B) <- colnames(P)

 colnames(A) <- colnames(B) <- colnames(pk) <- paste('LB',1:K,sep='')

 val_func <- xab$value

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

 #  names(results) <- c('Composition data matrix',
 #                      'Expected budget',
 #                      'Residual matrix',
 #                      'Mixing parameters',
 #                      'Latent budgets',
 #                      'Budget proportions',
 #                      'Value of the -loglik function',
 #                      'Number of iteractions') 
 # 
 class(results) <- c("lba.mle.fe",
                     "lba.mle")

 invisible(results)

}  

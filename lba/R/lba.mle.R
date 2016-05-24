lba.mle <- function(obj     , 
                    A       , 
                    B       , 
                    K       , 
                    tolG    , 
                    tolA    , 
                    tolB    ,
                    itmax.unide,
                    itmax.ide,
                    trace.lba,
                    toltype , 
                    what,
                    ...) 
{
 # The default toltype is 'all'
 I <- nrow(obj)
 J <- ncol(obj)

 P <- obj/rowSums(obj)
 #
 # ------------------------------------------------------------------------------
 # creating random generated values for alpha(k|i)
 if (is.null(A)) {
  A <- matrix(runif(I * K), 
              ncol = K)
  A <- A/rowSums(A)
 }
 #
 # ------------------------------------------------------------------------------
 # creating random generated values for beta(j|k)
 if (is.null(B)) {
  B <- matrix(runif(J * K), 
              ncol = K)
  B <- t(t(B)/colSums(B))
 }

 iter_unide <- 0  # counter of interactions numbers

 mle_func <- function(obj,
                      I,
                      J,
                      K,
                      A,
                      B){

  pip <- rowSums(obj)/sum(obj) #this is pi+
  pjki <- rep(0,I*J*K) #this will become pijk/pi+ see van der Ark page 80
  m <- 0
  # this makes i=1,j=1,k=1; i=2,j=1,k=1; i=3,j=1,k=1...; i=2,j=1,k=1 and so on.
  for(k in 1:K) for(j in 1:J) for(i in 1:I) { 
   m <-  m +1
   pjki[m] <- A[i,k]*B[j,k] 
  }

  mi <- matrix(pjki, I)                                          
  pijk <- mi*pip #this is pijk see van der Ark page 80 
  nijk <- as.vector(pijk)*sum(obj) 
  val_function <- -sum(nijk * log(pjki)) # -loglikelihood function

 } 
 repeat {
  iter_unide <- iter_unide + 1

  if(trace.lba){
   cat('Unidentified iteration: ',iter_unide,'\n')
   cat('fval = ', signif(mle_func(obj,I,J,K,A,B),4), '\n') 
  }  
  #
  # -------------------------------------------------------------------------------
  # previous G2 to be used in |G2 - G2a|<tol
  qij <- A %*% t(B)
  kij <- obj/rowSums(obj)
  G2a <- 2 * sum(obj * log(kij/qij))
  akia <- A
  bjka <- B

  #
  # -------------------------------------------------------------------------------
  ini <- list()
  # ----------------------------------------------------------------------------
  # STARTING THE (expectation)E-STEP
  # -------------------------------------------------------------------------------
  for (k in 1:K) ini[[k]] <- c(A[, k], 
                               B[, k])
  #
  # ------------------------------------------------------------------------------
  # initial values list with vectors a1(I),b1(J);... ak(I),bk(J); ... aK(I),
  # bK(J) where a's are alfa values and b's are beta values. The lenghts of the
  # a's are I and of the b's are J. The elements of the list are c(a1,b1), ...
  # c(ak,bk),... c(aK,bK). 

  pijk <- lapply(ini, function(xx) outer(xx[1:I], 
                                         xx[(I + 1):(I + J)]) * (rowSums(obj)/sum(obj)))
  # The elements of the list pijk are matrices (I x J) for fixed k and the number
  # of elements is K. They are alpha(k|i)beta(j|k)*pi+ = pijk|i*pi+

  mpijk <- sapply(pijk, 
                  function(x) x, simplify = TRUE)
  # Transform pijk into a matrix I*JxK. The columns are the colmuns of each
  # element (matrix) one after the other, of the first then the second and so on.

  sij <- rowSums(mpijk)
  # The sum over k of alpha(k|i)beta(j|k), vector of lenght IxJ

  lnijk <- lapply(pijk, 
                  function(xx) obj * xx)
  # The product of N(ij)*pijk. It is a list of K matrices IxJ each

  nijk <- lapply(lnijk, 
                 function(x) x/sij)
  # The list containing the values of n(hat)(ijk). Each element is a matrix
  # IxJ and the number of elements of the list is K.

  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  #Starting computation for maximization step
  #-------------------------------------------------------------------------------        

  sjnijk <- sapply(nijk, 
                   function(x) rowSums(x), simplify = TRUE)
  # sum over j for i and k fixed. sjnijk is a IxK matrix

  sjknijk <- apply(sjnijk, 
                   1, 
                   function(x) sum(x))
  # sum over j and k for i fixed. sjknijk is a vector of lenght I

  sinijk <- sapply(nijk, 
                   function(x) colSums(x), simplify = TRUE)
  # sum over i for j and k fixed. sinijk is a JxK matrix

  sijnijk <- apply(sinijk, 
                   2, 
                   function(x) sum(x))
  # sum over i and j and k fixed. sjknijk is a vector of lenght K
  # -----------------------------------------------------------------------------
  # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # maximization step
  # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  A <- sjnijk/sjknijk
  # next values alphahat1 to be used in the next step of the iteration

  B <- t(t(sinijk)/sijnijk)
  # next values betahat1 to be used in the next step of the iteration

  # aki; bjk

  pij <- A %*% t(B)
  # Values of pihat(j|i) = the sum over k of alphahat(k|i)*betahat(j|k)


  piplus <- rowSums(obj)/sum(obj)
  # The margins of the unconditional proportions

  pk <- piplus %*% A
  # budget proportions

  G2 <- 2 * sum(obj * log(obj/(pij * rowSums(obj))))
  chi2 <- sum(((obj - pij * rowSums(obj))^2)/obj)
  # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  at <- max(abs(A - akia))
  bt <- max(abs(B - bjka))
  Gt <- abs(G2 - G2a)

  #
  # ------------------------------------------------------------------------------
  if(iter_unide > itmax.unide){
   warning("maximum number of the unidentified iteractions exceeded" )
   break 
  }

  if (toltype == "all" & Gt < tolG & at < tolA & bt < tolB) {
   break  #default
  } else 
   if (toltype == "G2" & Gt < tolG) {
    break
   } else 
    if (toltype == "ab" & at < tolA & bt < tolB) {
     break
    }

 }  #closing the repeat

 val_func <- mle_func(obj,I,J,K,A,B) # -loglikelihood function


 if(K==1){

  Aoi <- A

  Boi <- B

  iter_ide  <- 'not applicable'

 } else {

  AB_outerAB_inner <- identif(A,
                              B,
                              P,
                              K,
                              trace.lba,
                              itmax.ide,
                              what)

  Aoi <- AB_outerAB_inner[[1]]

  Boi <- AB_outerAB_inner[[2]]

  iter_ide <- AB_outerAB_inner[[3]]  

 }

 pij <- A %*% t(B)

 residual <- P - pij

 pip <- rowSums(obj)/sum(obj) #this is pi+ 

 pk <- pip %*% Aoi

 colnames(pk) <- colnames(A) <- colnames(B) <- paste('LB',1:K,sep='')

 rescB <- rescaleB(x=obj,
                   Aoi,
                   Boi)

 colnames(rescB) <- colnames(B)
 rownames(rescB) <- rownames(B) 

 res <- list(P,       
             pij,
             residual,
             A, 
             B,
             Aoi,
             Boi,
             rescB,
             pk,
             val_func,
             iter_unide,
             iter_ide)

 names(res) <- c('P',
                 'pij',
                 'residual',
                 'A',
                 'B',
                 'Aoi',
                 'Boi',
                 'rescB',
                 'pk',
                 'val_func',
                 'iter_unide',
                 'iter_ide')

 #  names(res) <- c('Composition data matrix',
 #                  'Expected budget',
 #                  'Residual matrix',
 #                  'Unidentified mixing parameters',
 #                  'Unidentified latent budgets',
 #                  ifelse(what=='outer',
 #                         'Outer extreme mixing parameters',
 #                         'Inner extreme mixing parameters'),
 #                  ifelse(what=='outer',
 #                         'Outer extreme latent budgets',
 #                         'Inner extreme latent budgets'),
 #                  'Rescaled latent budgets',
 #                  'Budget proportions',
 #                  'Value of the -loglik function',
 #                  'Number of unidentified iteractions',
 #                  'Number of identified iteractions')
 # 
 class(res) <- "lba.mle" 

 invisible(res)            
}



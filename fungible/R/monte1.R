#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Program author: Niels Waller
# updated April 19, 2011
#
# Compute multivariate nonnormal data by the method of:
# Vale, C. & Maurelli, V. (1983). Simulating multivariate 
# nonnormal distributions. Psychometrika, 48(3):465-471.
#
#  R <-matrix(.5,3,3); diag(R) <- 1
#  X <- monte1(seed=123, nvar=3, nsub=25000,cormat=R, skewvec=c(1,1,1),kurtvec=c(2,2,2) )
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




monte1<-function(seed,nvar, nsub,cormat,skewvec,kurtvec)
{

 set.seed(seed)  
 #   
        
#=================================================================#
# *****************Arguement definitions**************************#
 
#****************************************************************#
#****************************************************************#
#
#   Beware: Traveling beyond this point is dangerous!
#****************************************************************#
#
#*****************************************************************#
#
# DEFINE external functions: mkclusb, mkclusc, mkclusd
#
#
#*****************************************************************#
#
########################################################################
    call<-match.call()

 mkclusb.prg <- function(nsub, nvar, skewvec, kurtvec, seed, desired.cor, orientation, callnum)
    {
    #=====================================================================#
    # mkclusb.prg                                                         #
    #                                                                     #
    # This program calls external function mkclusc.prg and mkclusd.prg    # 
    #                                                                     #
    # a program for generating nonnormal multivariate data with           # 
    # specified correlation structure using equations described in        #
    # Vale, D. C. & Maurelli, V. A. (1983). Simulating multivariate       # 
    #    nonnormal distributions.  Psychometrika, 48, 465-471.            #
    #                                                                     #
    # nvar=number of variables                                            #
    # skewvec=vector of skewness values for nvar variables                #
    # kurtvec=vector of kurtosis values for nvar variables                #
    # seed = seed number for data generation                              #
    # desired.cor is a matrix containing the target correlation matrix    #
    # for nonnormal data                                                  #
    #                                                                     #
    # upon completion the function returns X (nsub x nvar) with non-normal#
    # data sampled from a population with correlation: desired.cor        # 
    #=====================================================================# 

    mkclusc.prg <- function(x, skew = 0, kurt = 0)
      {

      #======================================================================#
      # mkclusc.prg                                                          #
      # This function is minimized to determine the weights needed to        # 
      # transform the data to desired skewness and kurtosis                  #
      # See Vale and Maurelli (1983)                                         #
      #                                                                      #
      # b, c, and d are weights used to transform normal data                #
      # f, g, and h are the three nonlinear equations that must be solved    #
      #    to find b, c, and d                                               #
      # skew and kurt are the desired skewness and kurtosis values           #
      #======================================================================# 
       b <- x[1]
       c <- x[2]
       d <- x[3]   
       f <- (b^2 + 6 * b * d + 2 * c^2 + 15 * d^2 - 1)
       g <- 2 * c * (b^2 + 24 * b * d + 105 * d^2 + 2) - skew
       h <- 24 * (b * d + c^2 * (1 + b^2 + 28 * b * d) + d^2 * (12 + 48 * b * d + 141 * c^2 + 225 * d^2)) - kurt
           return((f^2 + g^2 + h^2))
      }
   #----------------------------------------------------------------------#
   mkclusd.prg <- function(p, r, matr)
     {
     f <- (p * (matr[1, 2] * matr[2, 2] + 3 * matr[1, 2] * matr[2, 4] + 3 * matr[1, 4] * matr[2, 2] + 9 * matr[1, 4] * matr[
          2, 4]) + p^2 * (2 * matr[1, 3] * matr[2, 3]) + p^3 * (6 * matr[1, 4] * matr[2, 4])) - r
         return(f^2)
     }


   #=====================================================================#
   #  Generating weights for nonnormal transformation by minimizing      #
   #  a set of nonlinear equations                                       #
   #                                                                     #
   #  mkclusc.prg is called here                                         #
   #=====================================================================#
   # bcdvec is a vector of the b,c, and d weights from Vale and Maurelli (1983) #
   bcdvec <- matrix(0, nrow = nvar, ncol = 3)  #
   for(i in 1:nvar) {
     bcdvec[i,] <- optim(par = c(1.0, .0, .0),fn = mkclusc.prg,method="L-BFGS-B", lower=-2,upper=2, skew = skewvec[i], 
                   kurt = kurtvec[i],,control=list(ndeps=rep(1e-7,3)))$par
     bcdvec[i,1]<-bcdvec[i,1]*-1
     bcdvec[i,3]<-bcdvec[i,3]*-1          
                          
   }
   # (avec=a) a = -c in the Vale and Maurelli equation          
   avec <-  - bcdvec[, 2]  #

   #matrix of weights (the regression constants) 
   #used for transformation; (a,b,c,d) is a nvar X 4 (abcd) matrix#
        constant.mat <- as.matrix(cbind(avec, bcdvec))  #

   #================================================================#
   # intermediate correlation matrix (rxx)for normal data generation#
   # needed to determine the required correlation structure that    #
   # will result in the desired correlation structure after the data#
   # have been transformed to desired skewness and kurtosis         #
   # This is done by minimizing the function found in mkclusd.prg   #
   #                                                                #
   # matr contains the a b c d weights for the two variables of     #
   # desired cor                                                    #
   #================================================================#
   rxx <- desired.cor
   if(orientation == TRUE) {
     rxx <- matrix(c(rep(0, (nvar * nvar))), nrow = nvar, ncol = nvar)
        for(r in 2:nvar) {
          for(col in 1:(r - 1)) {
             rxx[r, col] <- optim(method="L-BFGS-B", par = 0.1, fn = mkclusd.prg,control=list(ndeps=1e-7),lower=-.99,upper=.99,r = desired.cor[r, col], 
                  matr = matrix(c(constant.mat[r,  ], constant.mat[col,  ]), nrow = 2, ncol = 4, byrow = TRUE))$par
             rxx[col, r] <- rxx[r, col]
          }
       }
            I <- diag(nvar) #
      # I is an identity matrix (nvar X nvar) used to place ones on the diagonal
      # of rxx
    rxx <- rxx + I  #
  }
  #======================================================================#
  #  The intermediate correlation matrix is decomposed using a           #
  #  Cholesky decomposition.  The resulting triangular matrix and stand- #
  #  ardized random variates are used to create a matrix of              #
  #  normally distributed data with the desired pattern of               #
  #  correlations                                                        #
  #======================================================================#
  set.seed(seed)  # within.clus.dev is an nsub X nvar data matrix used to develop clustered data#
  within.clus.dev <- matrix(rnorm(nvar * nsub), ncol = nvar, nrow = nsub) #
  within.clus.dev <- apply(within.clus.dev, 2, scale)
  floadt <- (chol(rxx))   #
  # floadt is the triangular matrix of loadings resulting from a Cholesky
  # decomposition of the intermediate correlation matrix
  X <- within.clus.dev %*% floadt #
  # X is data which has desired correlation structure (based on rxx) 
  #===============================================================#
  # The weight matrix is now used to transform the data to have   #
  # desired skewness, kurtosis and correlation structure          #
  # using Fleischman and Vale & Maurelli's power method           #
  # Fleishman, A.I (1978). A method for simulating non-normal     #
  # distributions, Psychometrika, 43, 521-532.                    #
  #===============================================================#
  one <- rep(1, nsub) # a vector of ones #
  X2 <- X^2
  X3 <- X^3
   for(i in 1:nvar) {
     ytemp <- cbind(one, X[, i], X2[, i], X3[, i])
     w <- matrix(c(constant.mat[i,  ]), ncol = 1)
     X[, i] <- ytemp %*% w
  }
  X <- apply(X, 2, scale) #
  # X contains the nonnormal, standardized data #
  return(X)
}




#==============================================================#
# mkclusd.prg                                                  #
# function based on equation found in Vale and Maurelli(1983)  #
# this function is minimized to obtain an intermediate         #
# correlation matrix (used for normally distributed data)      #
# which will work with certain weights to achieve the desired  #
# correlation structure for the transformed (nonnormal) data   #
#                                                              #
# p is the correlation that is being solved for (the intermed- #
#   iate correlation).  Start values are provided.             #
# r is provided. It is the correlation desired for the non-    #
#   normal data.                                               #
# matr is a matrix of transformation weights (a,b,c,d; from    #
#   constant.mat in mkclusb.prg                                #
#==============================================================#
#*****************************************************************
#
#                        MAIN PROGRAM STARTS HERE  
      
    
    callnum <- 1    #used by mkclusb.prg to inform user of current cluster#
 
    data <- mkclusb.prg(nsub, nvar, skewvec, kurtvec,seed = seed, desired.cor = cormat, orientation = 
                    TRUE, callnum = callnum)
                    
    result<-list(data = data, call=call, nsub=nsub, nvar=nvar, cormat=cormat,
                 skewvec=skewvec,kurtvec=kurtvec, seed=seed)
    class(result)<-"monte1"
    result
                       
}

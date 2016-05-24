monte<-function(seed=123,nvar=4, nclus=3, clus.size=c(50,50,50), 
               eta2=c(0.619, 0.401, 0.941, 0.929), cor.list=NULL, random.cor=FALSE,
               skew.list=NULL ,kurt.list=NULL,secor=NULL,compactness=NULL,sortMeans=TRUE)
{

    set.seed(seed)  #
    if(is.null(compactness)) compactness<-rep(1,nclus)
    
    if(!is.null(cor.list) || random.cor==TRUE)
        orientation<-TRUE
    else
        orientation<-FALSE    
        
    if(!is.null(skew.list) && !is.null(kurt.list))
       transform<-TRUE
    else
       transform<-FALSE 
 
       
    seed.vec <- c(25,43,47,21,38,42,18,49,41,22,30,3,33,27,45,16,35,37,2,5,20,26,31,29,12,44,
                  11,36,46, 8,10,50,15,34,39, 4,48,19,13,23,9,28,6,17,14,24,7,1,40,32)

      
        
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
    #library(stats)
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
#
        #print(cat("\n", "Generating transformed data for cluster", callnum, "\n")) #   
#=====================================================================#
#  Generating weights for nonnormal transformation by minimizing      #
#  a set of nonlinear equations                                       #
#                                                                     #
#  mkclusc.prg is called here                                         #
#=====================================================================#
# bcdvec is a vector of the b,c, and d weights from Vale and Maurelli (1983) #

        bcdvec <- matrix(0, nrow = nvar, ncol = 3)  #
        for(i in 1:nvar) {
            bcdvec[i,] <-  optim(par = c(1.0, .0, .0),fn = mkclusc.prg,method="L-BFGS-B", lower=-2,upper=2, skew = skewvec[i], 
                kurt = kurtvec[i],,control=list(ndeps=rep(1e-7,3)))$par
                        bcdvec[i,1]<-bcdvec[i,1]*-1
                        bcdvec[i,3]<-bcdvec[i,3]*-1          
                          
        }
# (avec=a) a = -c in the Vale and Maurelli equation          
        avec <-  - bcdvec[, 2]  #
#matrix of weights (the regression constants) 
#used for transformation; (a,b,c,d) is a nvar X 4 (abcd) matrix#
        constant.mat <- as.matrix(cbind(avec, bcdvec))  #
#
#================================================================#
# intermediate correlation matrix (rxx)for normal data generation#
# needed to determine the required correlation structure that    #
# will result in the desired correlation structure after the data#
# been transformed to desired skewness and kurtosis              #
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
#
#
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
#=================================================================#
# use this transformation if random cluster volumes are desired
#       compactness <- runif(length(clus.size))
#       compactness <- sqrt(length(clus.size) * (compactness/sum(compactness)))   #
#=================================================================# 
   
    if(transform==TRUE){
        skewvec <- matrix(0, nclus, nvar)
        kurtvec <- matrix(0, nclus, nvar)
         for(i in 1:nclus) {
            skewvec[i,  ] <- skew.list[[i]]
            kurtvec[i,  ] <- kurt.list[[i]]
         }
   }      
       
    
    callnum <- 1    #used by mkclusb.prg to inform user of current cluster#
    cornum <- round(1/secor^2 + 3, 0)
    nsub <- sum(clus.size)    #number of subjects #   
    sum.size <- c(0, cumsum(clus.size))   #
#=================================================================#
# sum.size is a nclus+1 vector with element j equal to  
# the cummulative cluster sample sizes 
#
#================See Equation 15 in Waller, Underhill, & Kaiser=================#
    eta2 <- 1 + (nsub * eta2 - nsub)/(eta2 * sum(clus.size * compactness^2) - (nsub * eta2) + nsub)   #
#==================================================================
# 
#=======make design.mat===and dummy factor scores===============#   
    design.mat <- matrix(0, nrow = nclus, ncol = nclus - 1)
    for(i in 1:(nclus - 1)) {
        design.mat[i, i] <- 1
    }
    dm <- design.mat    #
#
# dm or design.mat is a nclus X (nclus-1) matrix with all elements#
# 0 except for those where column and row number are equal (these are 1) #
# 
# initiate factor score matrix #
    fs <- matrix(0, nrow = nsub, ncol = nclus - 1)
    k <- 1  # put in appropriate dummy scores #
    for(i in 1:nclus) {
        for(j in 1:clus.size[i]) {
            fs[k,  ] <- dm[i,  ]
            k <- k + 1
        }
    }
#
#
# fs is the "dummy" factor score matrix (nsub X nclus-1) with 0's and 1's #
# used to demarcate cluster membership#
#
    rawfs <- fs # save backup of dummy factor scores
#===================standardize fs and orthogonalize=========#
    if(ncol(fs) == 1) {
        fs <- scale(fs)
    }
    if(ncol(fs) >= 2) {
        fs <- apply(fs, 2, scale)   #
        fs <- apply(princomp(fs, cor = TRUE)$scores, 2, scale)
    }
#
#===generate dummy factor loading matrix=====================#
#
    fl <- matrix(runif(nvar * (nclus - 1), min = -1, max = 1), nrow = nvar, ncol = nclus - 1)   #
# fl is a factor loading matrix(nvar X nclus-1) from random uniform variates# 
#
    sigmahat <- fl %*% t(fl)    #
#
#=================================================================#
# sigmahat is the sum of squares and crossproducts of the loadings
# 
# Note, at this point the communalities can be greater than 1, thus
# the factor loading matrix must be normalized by rows
#=================================================================#
    h2 <- diag(sigmahat)    #
#
# h2 are the sum of squared factor loadings #
#===========normalize rows of dummy factor loading matrix=====#
    signfl <- sign(fl)  #
#
# signfl is a matrix containing +1 and -1 for sign of loadings#
#
# this normalizes rows #
    for(i in 1:nvar) {
        fl[i,  ] <- sqrt(eta2[i]) * sqrt(fl[i,  ]^2/h2[i])
    }
# putting the correct sign back on the altered  fl#
    fl <- signfl * fl   #
#   
#=================================================================#
    sigmahat <- fl %*% t(fl)    #
    h2 <- diag(sigmahat)    #
#
#==============lmn=latent means==================================#
# the latent means (cluster centroids) are computed by post multiplying
# the normalized factor loading matrix (fl) by the standardized (and 
# orthogonalized) dummy factor scores
#================================================================#
    lmn <- matrix(0, nrow = nvar, ncol = nclus)
    for(k in 1:nvar) {
        j <- 1
        for(i in 1:nclus) {
            lmn[k, i] <- sum(fl[k,  ] * fs[j,  ])
            j <- j + clus.size[i]
        }
    }
    lmn <- t(lmn)
    if(sortMeans) lmn <- apply(lmn, 2, sort)    #
#================================================================#
#===================For zero nuisance covariance=================#
    if(orientation == FALSE) {
        within.clus.dev <- matrix(rnorm(nsub * nvar), nrow = nsub, ncol = nvar) #
#
# within.clus.dev is an nsub X nvar matrix containing the simulated data #
# without the latent means, i.e., the within cluster deviation scores #
#
        for(i in 1:nclus) {
            if(transform == TRUE) {                   
                       
#==============nonlinear transformation of data=======================#
#              external prog. mkclusb.prg called here                 #
#              the desired correlation matrix is an identity matrix   #
#              when orientation = FALSE and thus no transformation of   # 
#              cors is necessary                                      #
#=====================================================================#
                desired.cor <- diag(nvar)
                within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- mkclusb.prg(nsub = clus.size[i], nvar = nvar, skewvec
                   = skewvec[i,  ], kurtvec = kurtvec[i,  ], seed = seed + i, desired.cor = desired.cor, orientation = 
                  orientation, callnum = callnum)
                callnum <- callnum + 1
            }
            if(nvar > 1) {
                within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- compactness[i] * apply(within.clus.dev[(1 + sum.size[
                  i]):sum.size[i + 1],  ], 2, scale) %*% diag(sqrt(1 - h2)) #
            }
            if(nvar == 1) {
                within.clus.dev[(1 + sum.size[i]):sum.size[i + 1]] <- compactness[i] * scale(within.clus.dev[(1 + sum.size[i]):
                  sum.size[i + 1]]) * sqrt(1 - h2)
            }
        }
    }
#========================for nuisance covariance =================#
# Using the Cholesky decomposition of the desired correlation     #
# matrix we generate a loading matrix (possibly for each cluster) #
# for the within cluster scores.                                  #
# This is done cluster by cluster below  with either randomly     #
# chosen within cluster correlation structure (random.cor=TRUE)      #
# or with desired correlations set by the user (cor.list,         #
# which is renamed desired.cor).                                  #
#=================================================================#
    if(orientation == TRUE) {
       
#
# uscor is a matrix of random normal deviates that will equal 
# within.clus.dev after appropriate weighting by uf (uf=loading matrix 
# from random cor structure
        uscor <- matrix(rnorm(nsub * nvar), nrow = nsub, ncol = nvar)
        for(i in 1:nclus) {               
            uscor[(1 + sum.size[i]):sum.size[i + 1],  ] <- apply(uscor[(1 + sum.size[i]):sum.size[i + 1],  ], 2, scale)
        }
#
#============For within group correlations that are chosen randomly with
#============a standard error = secor
        if(random.cor == TRUE) {
            detx <- -99 #used for determinant of x
            while(detx < 0.05) {
                x <- matrix(rnorm(cornum * nvar), nrow = cornum, ncol = nvar)
                detx <- prod(eigen(abs(cor(x)))$values)
            }
#
# a matrix of x values, standard error of the absolute value of 
# correlations of x=secor
            uf <- chol(abs(cor(x)))
            uf <- t(uf) #
# uf is a triangular matrix from the Cholesky decomp. of abs(correlation of x)#
# (the target correlation matrix) #
            print(cat("Expected within group correlation matrix", "\n"))
            expect.cor <- round(uf %*% t(uf), 3)
            print(expect.cor)
            within.clus.dev <- matrix(0, nrow = nvar, ncol = nsub)
            uscor <- t(uscor)
            for(i in 1:nclus) {
#
#=====================================================================#
                if(transform == TRUE) {                             
                  if(transform == TRUE & i == 1) {
                    within.clus.dev <- t(within.clus.dev)
                  }
#==============================nonlinear transformation of data=====#
#==============================external prog. mkclusb.prg called here=#
                  within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- mkclusb.prg(nsub = clus.size[i], nvar = nvar, skewvec
                     = skewvec[i,  ], kurtvec = kurtvec[i,  ], seed = seed + i, desired.cor = expect.cor, orientation = 
                    orientation, callnum = callnum)
                  callnum <- callnum + 1
                }
#
                if(transform != TRUE) {
                  within.clus.dev[, (1 + sum.size[i]):sum.size[i + 1]] <- uf %*% uscor[, (1 + sum.size[i]):sum.size[i + 1]]
                }
            }
        }
#
#============If within cluster correlations are specified by the user==========#
#============via  (cor.list)===============================================#
        if(random.cor == FALSE) {
            within.clus.dev <- matrix(0, nrow = nvar, ncol = nsub)
            uscor <- t(uscor)
            for(i in 1:nclus) {
                if(transform == TRUE) {
                  print(i)
                  if(transform == TRUE & i == 1) {
                    within.clus.dev <- t(within.clus.dev)
                  }
                  desired.cor <- matrix(unlist(cor.list[i]), nrow = nvar)   # unlisting cor.list
#==============================nonlinear transformation of data================#
#==============================external prog. mkclusb.prg called here==========#
                  within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- mkclusb.prg(nsub = clus.size[i], nvar = nvar, skewvec
                     = skewvec[i,  ], kurtvec = kurtvec[i,  ], seed = seed + seed.vec[i], desired.cor = desired.cor, orientation = 
                    orientation, callnum = callnum)
                  callnum <- callnum + 1
                }
#
                if(transform != TRUE) {
                  uf <- chol(matrix(unlist(cor.list[i]), nrow = nvar, byrow = T))
                  uf <- t(uf)
                  within.clus.dev[, (1 + sum.size[i]):sum.size[i + 1]] <- uf %*% uscor[, (1 + sum.size[i]):sum.size[i + 1]]
                }
            }
            if(transform == FALSE) {
                within.clus.dev <- t(within.clus.dev)
            }
        }
#
        for(i in 1:nclus) {
            if(orientation == TRUE & random.cor == TRUE & transform == FALSE & i == 1) {
                within.clus.dev <- t(within.clus.dev)
            }
            within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- apply(within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ], 
                2, scale)
        }
        within.clus.dev <- within.clus.dev %*% diag(sqrt(1 - h2))
        for(i in 1:nclus) {
            within.clus.dev[(1 + sum.size[i]):sum.size[i + 1],  ] <- compactness[i] * within.clus.dev[(1 + sum.size[i]):sum.size[i +
                1],  ]
        }
    }
#=====End "if orientation =TRUE loop=====================================#
# now add in the mean structure ======================================= #
    data.means <- fs %*% t(fl)  
    id <- rep(seq(1, nclus, 1), clus.size)
 
 
 ## SORT MEANS == TRUE
    if(sortMeans){
        data.means <- apply(data.means, 2, sort)   
        dm<-apply(apply(data.means,2,unique),2,sort)
        for(i in 1:nclus){
           for(j in 1:nvar){ 
              data.means[id==i,j]<-dm[i,j]
           }   
        }
    }    
 
              
    data <- data.means + within.clus.dev    #fs <- apply(fs, 2, unique)

    data <- cbind(id, data)
    result<-list(data = data, lmn = lmn, fl = round(fl, 4), fs = fs,call=call,nclus=nclus,nvar=nvar,cor.list=cor.list,
                 skew.list=skew.list,kurt.list=kurt.list,clus.size=clus.size, eta2=eta2,seed=seed)
    class(result)<-"monte"
    result
}

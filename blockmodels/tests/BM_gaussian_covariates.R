require('blockmodels')
##
## SBM
##

## generation of one SBM network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
Y1 <- matrix(runif(n*n),n,n)
Y2 <- matrix(runif(n*n),n,n)
M<-matrix(rnorm(n*n,sd=5),n,n)+Z%*%Mu%*%t(Z)+4.2*Y1-1.6*Y2 ## adjacency matrix

## estimation
my_model <- BM_gaussian_covariates("SBM",M,list(Y1,Y2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
Mu[lower.tri(Mu)]<-t(Mu)[lower.tri(Mu)]
Y1 <- matrix(runif(n*n),n,n)
Y2 <- matrix(runif(n*n),n,n)
Y1[lower.tri(Y1)]<-t(Y1)[lower.tri(Y1)]
Y2[lower.tri(Y2)]<-t(Y2)[lower.tri(Y2)]
M<-matrix(rnorm(n*n,sd=5),n,n)+Z%*%Mu%*%t(Z)+4.2*Y1-1.6*Y2 ## adjacency matrix
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_gaussian_covariates("SBM_sym",M,list(Y1,Y2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

##
## LBM
##

## generation of one LBM network
npc <- c(20,10) # nodes per class
Q <- c(1,2) # classes
n <- npc * Q # nodes
Z1<-diag(Q[1])%x%matrix(1,npc[1],1)
Z2<-diag(Q[2])%x%matrix(1,npc[2],1)
Mu<-20*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
Y1 <- matrix(runif(n[1]*n[2]),n[1],n[2])
Y2 <- matrix(runif(n[1]*n[2]),n[1],n[2])
M<-matrix(rnorm(n[1]*n[2],sd=5),n[1],n[2])+Z1%*%Mu%*%t(Z2)+4.2*Y1-1.6*Y2 ## adjacency matrix

## estimation
my_model <- BM_gaussian_covariates("LBM",M,list(Y1,Y2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

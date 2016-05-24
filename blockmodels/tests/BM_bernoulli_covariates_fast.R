require('blockmodels')
##
## SBM
##

## generation of one SBM network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
sigmo <- function(x){1/(1+exp(-x))}
Z<-diag(Q)%x%matrix(1,npc,1)
Mg<-8*matrix(runif(Q*Q),Q,Q)-4
Y1 <- matrix(runif(n*n),n,n)-.5
Y2 <- matrix(runif(n*n),n,n)-.5
M_in_expectation<-sigmo(Z%*%Mg%*%t(Z) + 5*Y1-3*Y2)
M<-1*(matrix(runif(n*n),n,n)<M_in_expectation)

## estimation
my_model <- BM_bernoulli_covariates_fast("SBM",M,list(Y1,Y2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)


##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
sigmo <- function(x){1/(1+exp(-x))}
Z<-diag(Q)%x%matrix(1,npc,1)
Mg<-8*matrix(runif(Q*Q),Q,Q)-4
Mg[lower.tri(Mg)]<-t(Mg)[lower.tri(Mg)]
Y1 <- matrix(runif(n*n),n,n)-.5
Y2 <- matrix(runif(n*n),n,n)-.5
Y1[lower.tri(Y1)]<-t(Y1)[lower.tri(Y1)]
Y2[lower.tri(Y2)]<-t(Y2)[lower.tri(Y2)]
M_in_expectation<-sigmo(Z%*%Mg%*%t(Z) + 5*Y1-3*Y2)
M<-1*(matrix(runif(n*n),n,n)<M_in_expectation)
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_bernoulli_covariates_fast("SBM_sym",M,list(Y1,Y2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)



##
## LBM
##

## generation of one LBM network
npc <- c(20,10) # nodes per class
Q <- c(1,2) # classes
n <- npc * Q # nodes
sigmo <- function(x){1/(1+exp(-x))}
Z1<-diag(Q[1])%x%matrix(1,npc[1],1)
Z2<-diag(Q[2])%x%matrix(1,npc[2],1)
Mg<-8*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])-4
Y1 <- matrix(runif(n[1]*n[2]),n[1],n[2])-.5
Y2 <- matrix(runif(n[1]*n[2]),n[1],n[2])-.5
M_in_expectation<-sigmo(Z1%*%Mg%*%t(Z2) + 5*Y1-3*Y2)
M<-1*(matrix(runif(n[1]*n[2]),n[1],n[2])<M_in_expectation)

## estimation
my_model <- BM_bernoulli_covariates_fast("LBM",M,list(Y1,Y2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

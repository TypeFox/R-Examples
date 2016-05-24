require('blockmodels')
#
# SBM
#

## generation of one SBM network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
L<-70*matrix(runif(Q*Q),Q,Q)
M_in_expectation<-Z%*%L%*%t(Z)
M<-matrix(
    rpois(
        length(as.vector(M_in_expectation)),
        as.vector(M_in_expectation))
    ,n,n)

## estimation
my_model <- BM_poisson("SBM",M , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
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
L<-70*matrix(runif(Q*Q),Q,Q)
L[lower.tri(L)]<-t(L)[lower.tri(L)]
M_in_expectation<-Z%*%L%*%t(Z)
M<-matrix(
    rpois(
        length(as.vector(M_in_expectation)),
        as.vector(M_in_expectation))
    ,n,n)
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_poisson("SBM_sym",M , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
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
L<-70*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M_in_expectation<-Z1%*%L%*%t(Z2)
M<-matrix(
    rpois(
        length(as.vector(M_in_expectation)),
        as.vector(M_in_expectation))
    ,n[1],n[2])

## estimation
my_model <- BM_poisson("LBM",M , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

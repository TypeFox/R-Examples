require('blockmodels')
##
## SBM
##

## generation of one SBM network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
P00<-matrix(runif(Q*Q),Q,Q)
P10<-matrix(runif(Q*Q),Q,Q)
P01<-matrix(runif(Q*Q),Q,Q)
P11<-matrix(runif(Q*Q),Q,Q)
SumP<-P00+P10+P01+P11
P00<-P00/SumP
P01<-P01/SumP
P10<-P10/SumP
P11<-P11/SumP
MU<-matrix(runif(n*n),n,n)
M1<-1*(MU>Z%*%(P00+P01)%*%t(Z))
M2<-1*((MU>Z%*%P00%*%t(Z)) & (MU<Z%*%(P00+P01+P11)%*%t(Z))) ## adjacency matrices


## estimation
my_model <- BM_bernoulli_multiplex("SBM",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
P00<-matrix(runif(Q*Q),Q,Q)
P10<-matrix(runif(Q*Q),Q,Q)
P01<-matrix(runif(Q*Q),Q,Q)
P11<-matrix(runif(Q*Q),Q,Q)
SumP<-P00+P10+P01+P11
P00<-P00/SumP
P01<-P01/SumP
P10<-P10/SumP
P11<-P11/SumP
P00[lower.tri(P00)]<-t(P00)[lower.tri(P00)]
P01[lower.tri(P01)]<-t(P01)[lower.tri(P01)]
P10[lower.tri(P10)]<-t(P10)[lower.tri(P10)]
P11[lower.tri(P11)]<-t(P11)[lower.tri(P11)]
MU<-matrix(runif(n*n),n,n)
MU[lower.tri(MU)]<-t(MU)[lower.tri(MU)]
M1<-1*(MU>Z%*%(P00+P01)%*%t(Z))
M2<-1*((MU>Z%*%P00%*%t(Z)) & (MU<Z%*%(P00+P01+P11)%*%t(Z))) ## adjacency matrices


## estimation
my_model <- BM_bernoulli_multiplex("SBM_sym",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
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
P00<-matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
P10<-matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
P01<-matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
P11<-matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
SumP<-P00+P10+P01+P11
P00<-P00/SumP
P01<-P01/SumP
P10<-P10/SumP
P11<-P11/SumP
MU<-matrix(runif(n[1]*n[2]),n[1],n[2])
M1<-1*(MU>Z1%*%(P00+P01)%*%t(Z2))
M2<-1*((MU>Z1%*%P00%*%t(Z2)) & (MU<Z1%*%(P00+P01+P11)%*%t(Z2))) ## adjacency matrices


## estimation
my_model <- BM_bernoulli_multiplex("LBM",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

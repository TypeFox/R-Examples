
## 3 nodes, A, B, C, with two levels each.
A <- factor(NA,levels=paste("A",1:2,sep=""))
B <- factor(NA,levels=paste("B",1:2,sep=""))
C <- factor(NA,levels=paste("C",1:2,sep=""))
sim.df <- data.frame(A,B,C)

## graph: A|B,C  C|B  B
nw <- network(sim.df,specifygraph=FALSE,doprob=FALSE)
nw <- getnetwork(insert(nw,2,1,nocalc=TRUE))
nw <- getnetwork(insert(nw,3,1,nocalc=TRUE))
nw <- getnetwork(insert(nw,2,3,nocalc=TRUE))

## setup a proposal simprob, and correct it
sim.nw <- makesimprob(nw)
nodes(sim.nw)[[1]]$simprob[1:8] <- c(1/6,## P(A=1|B=1,C=1)=1/6
                                    5/6,## P(A=2|B=1,C=1)=5/6
                                    1/2,## P(A=1|B=2,C=1)=1/2
                                    1/2,## P(A=2|B=2,C=1)=1/2
                                    5/6,## P(A=1|B=1,C=2)=5/6
                                    1/6,## P(A=2|B=1,C=2)=1/6
                                    1/3,## P(A=1|B=2,C=2)=1/3
                                    2/3## P(A=2|B=2,C=2)=2/3
                                    )
nodes(sim.nw)[[2]]$simprob[1:2] <- c(1/2,## P(B=1)=1/2
                                    1/2## P(B=2)=1/2
                                    )

nodes(sim.nw)[[3]]$simprob[1:4] <- c(1/3,## P(C=1|B=1)=1/3
                                    2/3,## P(C=2|B=1)=2/3
                                    7/8,## P(C=1|B=2)=7/8
                                    1/8## P(C=2|B=2)=1/8
                                    )
## number of cases to simulate
n <- 1000

## do the simulation
set.seed(189)
sim <- rnetwork(sim.nw,n)

#### -- analysis 

## prior network 
simprior.nw <- network(sim)

## Imaginary database and joint parameter prior
N <- 500
simprior <- jointprior(simprior.nw,N)

res <- nwfsort(getnetwork(networkfamily(sim,simprior.nw,simprior)))

plot(res)

## remove equivalent networks
plot( unique(res,equi=TRUE) )

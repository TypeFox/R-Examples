
library(geigen)
source("testgsvd.R")

# from NAG documentation F08VNF
# complex
A <- matrix( c( .96-.81i , -.03+.96i , -.91+2.06i, -.05+.41i,
               -.98+1.98i, -1.20+.19i, -.66+.42i , -.81+.56i,
                .62-.46i ,  1.01+.02i,  .63-.17i ,-1.11+.6i,
                .37+.38i ,  .19-.54i ,-.98-.36i  ,  .22-.2i,
                .83+.51i ,  .2+.01i  ,-.17-.46i  , 1.47+1.59i,
               1.08-.28i ,  .2-.12i  ,-.07+1.23i ,  .26+.26i),
               nrow=6, ncol=4, byrow=TRUE)

B <- matrix(c(1,0,-1,0,0,1,0,-1), ncol=4, byrow=TRUE)

z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(18),2), nrow=3, byrow=TRUE)
B <- matrix(round(runif(30),2), ncol=6, byrow=TRUE)
A <- A+0i
B <- B+0i
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(18),2), nrow=3, byrow=TRUE)
B <- matrix(round(runif(24),2), ncol=6, byrow=TRUE)
A <- A+0i
B <- B+0i
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(36),2), ncol=6, byrow=TRUE)
B <- matrix(round(runif(36),2), ncol=6, byrow=TRUE)
A <- A+0i
B <- B+0i
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(36),2), ncol=6, byrow=TRUE)
B <- matrix(round(runif(24),2), ncol=6, byrow=TRUE)
A <- A+0i
B <- B+0i
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
A <- matrix(round(runif(12),2), nrow=2, byrow=TRUE)
B <- matrix(round(runif(24),2), ncol=6, byrow=TRUE)
A <- A+0i
B <- B+0i
z <- gsvd(A,B)
testgsvd(z,A,B)

set.seed(11)
NrowA <- 25
NrowB <- 5
Ncol  <- 6
A <- matrix(round(runif(NrowA*Ncol),2), ncol=Ncol, byrow=TRUE)
B <- matrix(round(runif(NrowB*Ncol),2), ncol=Ncol, byrow=TRUE)
A <- A+matrix(round(runif(NrowA*Ncol),2)*1i, ncol=Ncol, byrow=TRUE)
B <- B+matrix(round(runif(NrowB*Ncol),2)*1i, ncol=Ncol, byrow=TRUE)
z <- gsvd(A,B)
testgsvd(z,A,B)


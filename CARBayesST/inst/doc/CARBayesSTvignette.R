### R code from vignette source 'CARBayesSTvignette.Rnw'

###################################################
### code chunk number 1: CARBayesSTvignette.Rnw:262-263
###################################################
library(CARBayesST)


###################################################
### code chunk number 2: CARBayesSTvignette.Rnw:329-335
###################################################
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)
N <- 10
N.all <- N * K


###################################################
### code chunk number 3: CARBayesSTvignette.Rnw:340-349
###################################################
W <-array(0, c(K,K))
    for(i in 1:K)
    {
        for(j in 1:K)
        {
        temp <- (Grid[i,1] - Grid[j,1])^2 + (Grid[i,2] - Grid[j,2])^2
            if(temp==1)  W[i,j] <- 1 
        }    
    }


###################################################
### code chunk number 4: CARBayesSTvignette.Rnw:354-362
###################################################
D <-array(0, c(N,N))
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            if(abs((i-j))==1)  D[i,j] <- 1 
        }    
    }


###################################################
### code chunk number 5: CARBayesSTvignette.Rnw:367-368
###################################################
Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))


###################################################
### code chunk number 6: CARBayesSTvignette.Rnw:373-375
###################################################
Q.W.inv <- solve(Q.W)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.01 * Q.W.inv))


###################################################
### code chunk number 7: CARBayesSTvignette.Rnw:381-384
###################################################
Q.D <- 0.99 * (diag(apply(D, 2, sum)) - D) + 0.01 * diag(rep(1,N))
Q.D.inv <- solve(Q.D)
delta <- mvrnorm(n=1, mu=rep(0,N), Sigma=(0.01 * Q.D.inv))


###################################################
### code chunk number 8: CARBayesSTvignette.Rnw:389-394
###################################################
phi.long <- rep(phi, N)
delta.long <- kronecker(delta, rep(1,K))
LP <- 4 + phi.long +  delta.long
mean <- exp(LP)
Y <- rpois(n=N.all, lambda=mean)



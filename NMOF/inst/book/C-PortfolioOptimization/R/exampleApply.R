# exampleApply.R -- version 2010-10-09
# objective function
omega <- function(r,theta) {
    rr <- r - theta
    omega <- -sum(rr - abs(rr)) / sum(rr + abs(rr))
    omega
}

# create artifical data 
# ...ns = number of scenarios 
# ...na = number of assets
ns <- 200L
na <- 100L
R  <- array(rnorm(ns*na)*0.05, dim = c(ns,na) )

# set up a random portfolio
w   <- runif(na)
w   <- w / sum(w)

# compute returns
rp <- R %*% w

# compute omega
omega(rp, theta = 0.001)

# objective function, alternative
omega2 <- function(r,theta) {
    rr <- r - theta
    omega2 <- -colSums(rr - abs(rr)) / colSums(rr + abs(rr))
    omega2
}

# check: compute omega
omega2(rp, theta = 0.001)

# set up a random population
# ...nP = population size
nP <- 200L
P  <- array(runif(na*nP), dim = c(na,nP))
P  <- P / outer(numeric(na) + 1, colSums(P))  # budget constraint

# evaluate population 
# ...variant 1
rp <- R %*% P
system.time({
    for (r in 1L:100L){
        for(i in 1L:nP) a1 <- omega(rp[ ,i], theta = 0.001) 
    }
})

# ...variant 2
rp <- R %*% P
system.time({
    for (r in 1L:100L) 
        a2 <- apply(rp, 2, omega, theta = 0.001)
})

# ...variant 3
rp <- R %*% P
system.time({
    for (r in 1L:100L) 
        a3 <- omega2(rp, theta = 0.001)
})
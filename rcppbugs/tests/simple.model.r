NR <- 100
NC <- 2
b <- normal(rnorm(NC),mu=0,tau=0.001)
X <- matrix(rnorm(NR * NC),NR,NC)
tau.y <- uniform(rnorm(1),0,100)

y.hat <- function() {
    X %*% b
}

y <- normal(rnorm(NR),observed=TRUE,mu=y.hat,tau=tau.y)

run.model(list(b,tau.y,y.hat,y), iterations=1e4, burn=1e3, adapt=1e3, thin=10)

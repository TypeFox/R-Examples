x <- 0:2
sum( (x - 1)^2 * dbinom(x,2,0.5) )   # same as above
n <- 5; p <- 0.2; x <- 0:n
sum( (x - n*p)^2 * dbinom(x,n,p) )   # X ~ Binom(5,0.2)
n <- 5; p <- 0.8; x <- 0:n
sum( (x - n*p)^2 * dbinom(x,n,p) )   # X ~ Binom(5,0.8)
n <- 10; p <- 0.8; x <- 0:n
sum( (x - n*p)^2 * dbinom(x,n,p) )   # X ~ Binom(10,0.8)
n <- 20; p <- 0.8; x <- 0:n
sum( (x - n*p)^2 * dbinom(x,n,p) )   # X ~ Binom(20,0.8)

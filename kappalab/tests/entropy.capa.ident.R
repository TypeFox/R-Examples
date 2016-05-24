library(kappalab)

## a set of random generated profiles
## for instance, marks on a [0,20] scale
p <- data.frame(matrix(runif(500,0,20),100,5))
names(p) <- c("Stat","Prob","Alg","Cal","Eng")

## discretization
p[p <= 5] <- 1  
p[p > 5 & p <= 10] <- 2 
p[p > 10 & p <= 15] <- 3 
p[p > 15] <- 4

d <- data.frame(factor(p[[1]]),
                factor(p[[2]]),
                factor(p[[3]]),
                factor(p[[4]]),
                factor(p[[5]]))

## associated unsupervised capacity
mu <- entropy.capa.ident(d)
mu




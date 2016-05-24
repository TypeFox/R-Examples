rm(list=ls())
library(quadrupen)
p <- 20
n <- 10
C <- bandSparse(p,k=0:1,diagonals=list(rep(1,p),rep(-1,p-1)))
out <- bounded.reg(matrix(rnorm(20*10),10,20),rnorm(10), struct=t(C) %*% C, lambda2=1)
plot(out)

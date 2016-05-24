# Test OPTIONS input (Unit test)
# This test example is provided by Hongbo Dong.

library(Rdsdp)

K=c()
K$l=0
K$s=c(3)
C = matrix(c(1,0,0,0,1,0,0,0,1),1,9, byrow=TRUE)
A = matrix(c(1,0,0,0,0,0,0,0,0),1,9,byrow=TRUE)
b <- c(1);

OPTIONS=c()
OPTIONS$print=1
OPTIONS$logsummary=1
OPTIONS$outputstats=1

result = dsdp(A,b,C,K,OPTIONS)

stopifnot(all.equal(result$y,c(1),tolerance=1e-05))

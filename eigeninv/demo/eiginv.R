n <- 100
evals <- sort(rnorm(n)) 
system.time(A <- eiginv(evals))
all.equal(evals, eigen(A)$val)

ev <- -seq(1:6)
set.seed(123)
B <- eiginv(ev)
eigen(B)$val
B2 <- eiginv(ev) # B and B2 will be different
all.equal(B, B2)

set.seed(123)
B3 <- eiginv(ev)
all.equal(B, B3) # will be identical

n <- 9
evals <- c(1, 0.95, rev(sort(runif(n-2, 0, 0.9))))
B4 <- eiginv(evals, stoch=TRUE)
eigen(B4)$value
rowSums(B4)

B5 <- eiginv(evals, symm=TRUE, stoch=TRUE)
eigen(B5)$value
rowSums(B5)
colSums(B5)


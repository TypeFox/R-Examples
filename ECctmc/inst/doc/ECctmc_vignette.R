## ----example1------------------------------------------------------------
set.seed(183427)
require(ECctmc)
# rates
r1 <- 1 # 1->2
r2 <- 0.75 # 2->3
r3 <- 0.5 # 3->1
r4 <- 0.5 # 3-> 2
Q <- matrix(c(-r1, r1, 0, 0, -r2, r2, r3, r4, -(r3+r4)), nrow = 3, byrow = TRUE)

# sample path
path <- sample_path(a=1, b=3, t0=0, t1=5, Q=Q)
path
plot(stepfun(x=path[1:(nrow(path)-1),"time"], y = path[,"state"]), xlim = c(0,5), xlab = "Time", ylab = "State", main = "Sample path")


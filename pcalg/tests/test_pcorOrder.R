library(pcalg)
## require(RBGL)

set.seed(123)

g <- randomDAG(4,0.6)
dat <- rmvDAG(1000, g, errDist="normal")

res1 <- pcorOrder(3,4,c(2,1), cor(dat))
res2 <- pcorOrder(3,4,  2,    cor(dat))

stopifnot(all.equal(res1, -0.0047417463),
          all.equal(res2,  0.0089093413))

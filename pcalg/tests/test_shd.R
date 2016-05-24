library(pcalg)

set.seed(123)

## expl 1
g1 <- randomDAG(5,0.5)
g2 <- randomDAG(5,0.5)
res1 <- (shd(g1,g2)==4)

## expl 2
g3 <- dag2cpdag(g1)
res2 <- (shd(g3,g1)==3)

if(!all(c(res1,res2))) {
  stop("Test of shd: Theoretical value not matched!")
}

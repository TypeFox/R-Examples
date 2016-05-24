library(EMCluster, quiet = TRUE)
set.seed(1234)

x <- da1$da
TC <- da1$class

n <- nrow(x)
p <- ncol(x)
k <- 10

ret.em <- init.EM(x, nclass = k, method = "em.EM")
ret.Rnd <- init.EM(x, nclass = k, method = "Rnd.EM", EMC = .EMC.Rnd)
ret.Rndp <- init.EM(x, nclass = k, method = "Rnd.EM", EMC = .EMC.Rndp)
ret.svd <- emgroup(x, nclass = k)

par(mfrow = c(2, 2))
plotem(ret.em, x, main = "em")
plotem(ret.Rnd, x, main = "Rnd")
plotem(ret.Rndp, x, main = "Rnd+")
plotem(ret.svd, x, main = "svd")

ret.all <-
cbind(
  c(ret.em$llhdval, ret.Rnd$llhdval, ret.Rndp$llhdval,
    ret.svd$llhdval),
  c(RRand(ret.em$class, TC)$adjRand,
    RRand(ret.Rnd$class, TC)$adjRand,
    RRand(ret.Rndp$class, TC)$adjRand,
    RRand(ret.svd$class, TC)$adjRand)
)
rownames(ret.all) <- c("em", "Rnd", "Rnd+", "svd")
colnames(ret.all) <- c("logL", "adjR")
ret.all

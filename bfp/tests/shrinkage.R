## test the shrinkage sampling.

library(bfp)

n <- 1e+5
R2 <- 0.7
nObs <- 100
p <- 20
alpha <- 3.5

tVec <- bfp:::rshrinkage(n=n,
                                    R2=R2,
                                    nObs=nObs,
                                    p=p,
                                    alpha=alpha)
histRet <- hist(tVec,
                nclass=50,
                prob=TRUE)

## compare with unn density:
grid <- histRet$mids

vals <- (1 - grid)^((p + alpha - 2) / 2 - 1) * (1 - R2 * grid)^(-(nObs - 1) / 2)

## scale
vals <- vals / max(vals) * max(histRet$density)

lines(grid,
      vals,
      col=2)

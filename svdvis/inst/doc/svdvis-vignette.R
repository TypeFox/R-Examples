## ---- fig.show='hold'----------------------------------------------------
set.seed(1234)
library(svdvis)
B = c(runif(100, min=0, max=1), rep(0,400))
L = c(rep(1, 10), rep(-1, 10))
L = L / sd(L)
E = matrix(rnorm(500*20), nrow=500)
Y = B %*% t(L) + E

svd.obj = svd(Y)
colnames(svd.obj$v) = paste0("V",1:20)
rownames(svd.obj$v) = paste0("Sample",1:20)

## ------------------------------------------------------------------------
svd.scree(svd.obj, subr=5,
          axis.title.x="Full scree plot", axis.title.y="% Var Explained")

## ------------------------------------------------------------------------
svd.scatter(svd.obj)

## ------------------------------------------------------------------------
svd.scatter(svd.obj, r=3, alpha=.5,
            group=c(rep("Group 1", 10), rep("Group 2", 10)))

## ------------------------------------------------------------------------
svd.heatmap(svd.obj, r=5)

## ------------------------------------------------------------------------
svd.parallel(svd.obj, r=5, alpha=.5,
             group=c(rep("Group 1", 10), rep("Group 2", 10)))

## ------------------------------------------------------------------------
svd.radial(svd.obj, r=3,
           group=c(rep("Group 1", 10), rep("Group 2", 10)))


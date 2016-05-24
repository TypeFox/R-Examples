## Check the equivalence between manual and automatic scaling
library(SOMbrero)
data(lesmis)

##
set.seed(321)
nsom1 <- trainSOM(iris[,1:4], type= "numeric", maxit= 10, scaling= "unitvar")
iris.scale <- scale(iris[,1:4], center= TRUE, scale= TRUE)
rownames(iris.scale) <- rownames(iris)
set.seed(321)
nsom2 <- trainSOM(iris.scale, type= "numeric", maxit= 10, scaling= "none")
stopifnot(identical(nsom1$clustering, nsom2$clustering))

##
set.seed(321)
nsom1 <- trainSOM(iris[,1:4], type= "numeric", maxit= 10, scaling= "center")
iris.scale <- scale(iris[,1:4], center= TRUE, scale= FALSE)
rownames(iris.scale) <- rownames(iris)
set.seed(321)
nsom2 <- trainSOM(iris.scale, type= "numeric", maxit= 10, scaling= "none")
stopifnot(identical(nsom1$clustering, nsom2$clustering))

## 
set.seed(321)
rsom1 <- trainSOM(dissim.lesmis, type= "relational", maxit= 10, 
                  scaling= "cosine", nb.save=10)
sim.matrix <- -.5* (diag(1, nrow(dissim.lesmis))-1/nrow(dissim.lesmis)) %*%
  dissim.lesmis %*% (diag(1, nrow(dissim.lesmis))-1/nrow(dissim.lesmis))
sim.matrix <- round(sim.matrix, 10)
scaled.ker <- sweep(sweep(sim.matrix,2,sqrt(diag(sim.matrix)),"/"),
                    1,sqrt(diag(sim.matrix)),"/")
scaled.diss <- 2-2*scaled.ker
colnames(scaled.diss) <- colnames(dissim.lesmis)
rownames(scaled.diss) <- rownames(dissim.lesmis)
lesmis.scale <- round(scaled.diss,10)
set.seed(321)
rsom2 <- trainSOM(lesmis.scale, type= "relational", maxit= 10, scaling= "none",
                  nb.save=10)
stopifnot(identical(rsom1$clustering, rsom2$clustering))


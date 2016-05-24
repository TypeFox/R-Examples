library("flexclust")

x <- matrix(runif(1000, 0, 1), ncol=2)

set.seed(1)
cl0 <- kcca(x, k=20)
cl1 <- kcca(x, k=20, weights=x[,1])
cl2 <- kcca(x, k=20, weights=x[,2])

pdf("weights.pdf")
plot(cl0)
plot(cl1)
plot(cl2)

m1 <- colMeans(cl1@centers)
m2 <- colMeans(cl2@centers)

stopifnot(diff(m1)<0)
stopifnot(diff(m2)>0)


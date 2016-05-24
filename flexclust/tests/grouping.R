library(flexclust)

set.seed(1207)

x <- cbind(sort(runif(1000)), runif(1000))

fam1 <- kccaFamily("kmeans")
fam2 <- kccaFamily("kmedians")
fam3 <- kccaFamily("kmedians", groupFun="differentClusters")

k1 <- 30
k2 <- 5

cl1 <- kcca(x, k=k1, family=fam1)
cl1
cl2 <- kcca(x, k=k2, group=predict(cl1), family=fam2)
cl2

tab12 <- table(predict(cl1),predict(cl2))
tab12
stopifnot(sum(tab12>0)==k1)

###**********************************************************

## the first column of x is sorted, hence recycling group membership
## gives observations in each group that are easy to seperate

g <- 1:(nrow(x)/10)
g <- rep(g, 10)
cl3 <- kcca(x, k=k1, group=g, family=fam3)
tab3 <- table(predict(cl3),g)
tab3
stopifnot(all(tab3 %in% 0:1))

###**********************************************************

cl4 <- stepFlexclust(as.matrix(iris[,-5]), k=10, group = iris[, 5])
tab4 <- table(predict(cl4),iris$Species)
tab4
stopifnot(all(tab4 %in% c(0, 50)))




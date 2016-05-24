data(disc2d)

### no partition and non parallel hmac ###
iris.hmac.nonparallel=phmac(iris[,-5],npart=1,parallel=FALSE)
### plot the hieararchical tree
plot(iris.hmac.nonparallel,level=3)

### parallel computing with partitioned data ###
disc2d.hmac.parallel=phmac(disc2d,npart=2,parallel=TRUE)

### View the summary of hmac object
summary(disc2d.hmac.parallel)

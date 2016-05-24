library("flexclust")

## make sure bootFlexclust works also for only one number of clusters
data(Nclus)
b <- bootFlexclust(Nclus, 3, nboot=10, multicore=FALSE)
print(b)
summary(b)

pdf("boot.pdf")
plot(b)
densityplot(b)
dev.off()

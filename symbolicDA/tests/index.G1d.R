# Example 1
library("symbolicDA")
library(stats)
data("cars",package="symbolicDA")
x<-cars
d<-dist.SDA(x, type="U_2")
wynik<-hclust(d, method="ward", members=NULL)
clusters<-cutree(wynik, 4)
G1d<-index.G1d(d, clusters)
print(G1d)

# Example 2


data("cars",package="symbolicDA")
md <- dist.SDA(cars, type="U_3", gamma=0.5, power=2)
min_nc=2
max_nc=10
res <- array(0,c(max_nc-min_nc+1,2))
res[,1] <- min_nc:max_nc
clusters <- NULL
for (nc in min_nc:max_nc)
{
cl2 <- pam(md, nc, diss=TRUE)
res[nc-min_nc+1,2] <- G1d <- index.G1d(md,cl2$clustering)   
clusters <- rbind(clusters, cl2$clustering)
}
print(sprintf("max G1d for %f clusters=%.3f",max(res[,2]),(min_nc:max_nc)[which.max(res[,2])]))
print("clustering for max G1d")
print(clusters[which.max(res[,2]),])
write.table(res,file="G1d_res.csv",sep=";",dec=",",row.names=TRUE,col.names=FALSE)
plot(res, type="p", pch=0, xlab="Number of clusters", ylab="G1d", xaxt="n")
axis(1, c(min_nc:max_nc))

library(symbolicDA)
library(rgl)
library(clusterSim)
means <- matrix(c(0,0,0,
0,0,6,
0,6,0,
0,6,6,
6,0,0,
6,0,6,
6,6,0,
6,6,6),8,3,byrow=TRUE)
means<-means*1.5
means[5:8,1]<-means[5:8,1]-2
means[5:8,3]<-means[5:8,3]-2
cov <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)
t<-cluster.Gen(model=2, means=means, cov=cov, dataType="s", numObjects=10)
plot3dInterval(t$data, colors=rainbow(8)[t$clusters])
rgl.viewpoint(15,20,30)
rgl.snapshot("8_clusters_3d.jpg")

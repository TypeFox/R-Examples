
recluster.test.dist <- function (mat1,mat2,member,perm=1000,elev=2)
{
	test<-matrix(NA,perm,3)	
	res<-NULL
	dist1<-sum(dist(mat1[,1:2])^elev)/length(as.vector(dist(mat1[,1:2])))
	dist2<-sum(dist(mat2[,1:2])^elev)/length(as.vector(dist(mat2[,1:2])))
	res$ratio<-dist2/dist1
	dist3<-sum(dist(mat1[,1])^elev)/length(as.vector(dist(mat1[,1:2])))
	dist4<-sum(dist(mat2[,1])^elev)/length(as.vector(dist(mat2[,1:2])))
	res$ratioX<-dist4/dist3
	dist5<-sum(dist(mat1[,2])^elev)/length(as.vector(dist(mat1[,1:2])))
	dist6<-sum(dist(mat2[,2])^elev)/length(as.vector(dist(mat2[,1:2])))
	res$ratioY<-dist6/dist5
	for (i in 1:perm){
		memb2<-sample(member)
		mat3<-recluster.group.col(mat1,memb2)
		test[i,1]<-(sum(dist(mat3$aggr[,1:2])^elev)/length(as.vector(dist(mat3$aggr[,1:2]))))/dist1
		test[i,2]<-(sum(dist(mat3$aggr[,1])^elev)/length(as.vector(dist(mat3$aggr[,1:2]))))/dist3
		test[i,3]<-(sum(dist(mat3$aggr[,2])^elev)/length(as.vector(dist(mat3$aggr[,1:2]))))/dist5
		}
points(res$ratio,0,cex=1,col="red")
res$test<-mean(test[,1]>=res$ratio)
res$testX<-mean(test[,2]>=res$ratioX)
res$testY<-mean(test[,3]>=res$ratioY)
return (res)
}
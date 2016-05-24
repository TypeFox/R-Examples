recluster.group.col <- function (mat,member)
{
	res<-NULL
	member<-as.numeric(member)
	mat<-as.matrix(mat)
	aggr <- aggregate(mat, by = list(member), FUN = mean)
	rownames(aggr)<-aggr[,1]
	aggr[,1]<-NULL
	aggr<-as.matrix(aggr)
	all<-mat
	for (i in 1:nrow(mat)){
		pos<-member[i]
		all[i,]<-aggr[pos,]
	}
	res$aggr<-aggr
	res$all<-all
	return (res)
}

# "get clusters" wrappers 

agnes.wrap <- function(data.tree, clust.num=0)
{
	return( cutree(data.tree, clust.num) )
}

diana.wrap = agnes.wrap
mona.wrap = agnes.wrap # other wrapper !!!
hclust.wrap = agnes.wrap

kmeans.wrap <- function(data.cls, clust.num=0)
{
	return( data.cls$cluster )
}

pam.wrap <- function(data.cls, clust.num=0)
{
	return( data.cls )
}

clara.wrap <- function(data.cls, clust.num=0)
{
	return( data.cls$clustering )
}

# cluster methods wrappers

agnes.clust <- function(data, clust.num=0, method="empty", ... )
{
#cat("\nagnes.clust")
	return( agnes(data, method=method, ... ) )
}

diana.clust <- function(data, clust.num=0, method="empty", ... )
{
#cat("\ndiana.clust")
	return( diana(data, ... ) )
}

mona.clust <- function(data, clust.num=0, method="empty", ... )
{
#cat("\nmona.clust")
	return( mona(data, ... ) )
}

hclust.clust <- function(data, clust.num=0, method="empty", ... )
{
#cat("\nhclust.clust")
	return( hclust( dist(data), method=method, ... ) )
}

kmeans.clust <- function(data, clust.num=0, method="empty", ... )
{
#cat("\nkmeans.clust")
	return( kmeans(data, centers=clust.num, ... ) )
}

pam.clust <- function(data, clust.num=0, method="empty", ... )
{
#cat("\npam.clust")
	return( pam(data, clust.num, cluster.only = TRUE,  ... ) )
}

clara.clust <- function(data, clust.num=0, method="empty", ... )
{
#cat("\nclara.clust")
	return( clara(data, clust.num, ... ) )
}

# classification wrappers

knn.pred <- function( base.data, base.class, rest.data )
{
#cat("\nknn.pred")
	return( as.integer( knn(base.data, rest.data, base.class) ) )
}
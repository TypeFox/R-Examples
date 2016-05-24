msBP.postCluster <- function(y, weights)
{
n <- length(y)
vec.weights <- tree2vec(weights)
maxS <- weights$max.s
res<-.C("postCluster_C", s=as.integer(rep(0,n)), h=as.integer(rep(1,n)), y=as.double(y), as.double(vec.weights), as.integer(maxS), as.integer(n), as.integer(0), PACKAGE = "msBP")
cbind(res$s,res$h)
}

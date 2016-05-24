print.pdclust <-
function(x, ...)
{
	cat(paste("Permutation Distribution Clustering\n\n", 
	"Embedding dimension: ",x$m, 
	"\nTime delay	    : ",x$t,
	"\nNumber of objects: ",x$N,
	"\nClustering method: ",x$method,"\n",sep=""));

	invisible(x)	
}

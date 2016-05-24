str.pdclust <- function(object, ...)
{
#	class(object) <- "hclust"
#	dd <- as.dendrogram(object)
#	str(dd, ...)
	str(as.dendrogram(object,...))
}
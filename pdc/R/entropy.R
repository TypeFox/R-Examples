entropy <-
function(dist)
{
	dist <- dist[dist!=0]
	return(-sum(dist*log(dist))/log(length(dist)));
	
}

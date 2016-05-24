calculate.empirical.AUC <-
function(data, marker, status, tag.healthy, confidence.level) {
      marker.diseased = data[data[,status] != tag.healthy, marker]
    	n.diseased = length (marker.diseased)

    	marker.healthy = data[data[,status] == tag.healthy, marker]
    	n.healthy = length(marker.healthy)

    	# Function that counts the number of elements in a vector that take a value of zero:
    	count.zeros <- function(x)
    	{
        	length(x[x == 0])
    	}

    	# Function that counts the number of elements in a vector that take a value lower than zero:
    	count.neg <- function(x)
    	{
        	length(x[x < 0])
    	}

    	marker.diseasedmat <- matrix(rep(marker.diseased,n.healthy), nrow = n.healthy, byrow = T)
    	marker.healthymat <- matrix(rep(marker.healthy,n.diseased), nrow = n.healthy, byrow = F)
    	diffmat <- marker.healthymat-marker.diseasedmat
    	area <- (length(diffmat[diffmat < 0])+0.5*length(diffmat[diffmat == 0]))/(n.diseased*n.healthy)
    	
    	res <- list(AUC = area)
    	return(res)
}

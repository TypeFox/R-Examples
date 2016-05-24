
# This function accepts a "dist" or matrix of scores and
# returns an ordering, based on hierarchical clustering.
# If reorder is FALSE, the order returned by hclust is used, 
# otherwise clusters are ordered by placing the nearest end points
# adjacent to each other at a merge.
 
order.hclust <- 
function(merit,reorder=TRUE,...) {
    dis <- - merit
    if (is.matrix(dis)) 
    disd <- as.dist(dis)
    else {
	disd <- dis
	dis <- as.matrix(dis)}
    n <- nrow(dis)
    if (n <= 2)
    ord <- 1:n
    else {
	hc <- hclust(disd,...)
	if (reorder)
	hc <- reorder.hclust(hc,dis)
	ord <- hc$order}
    ord }


# This function accepts hc, the results of a hierarchical clustering
# and a "dist" or distance matrix. It returns a hierarchical clustering obtained by placing 
# the nearest end points adjacent to each other at each
#  merge of the hierarchical clustering

reorder.hclust <-
function(x,dis,...) {
    if (! is.matrix(dis)) dis <- as.matrix(dis)
    merges <- x$merge
    n <- nrow(merges)
    endpoints <- matrix(0,n,2)
    dir <- matrix(1,n,2)
    for (i in 1:n) {
	j <- merges[i,1]
	k <- merges[i,2]
	if ((j < 0) && (k < 0)) {
	    endpoints[i,1] <- -j
	    endpoints[i,2] <- -k}
	else if (j < 0) {
	    j <- -j
	    endpoints[i,1] <- j
	    if (dis[j,endpoints[k,1]] < dis[j,endpoints[k,2]])        
	    endpoints[i,2] <- endpoints[k,2]
	    else {
		endpoints[i,2] <- endpoints[k,1]
		dir[i,2] <- -1}}
	else if (k < 0) {
	    k <- -k
	    endpoints[i,2] <- k     
	    if (dis[k,endpoints[j,1]] < dis[k,endpoints[j,2]]){
		endpoints[i,1] <- endpoints[j,2]
		dir[i,1] <- -1 }
	    else {
		endpoints[i,1] <- endpoints[j,1]
	    }}
	else {
	    d11 <- dis[endpoints[j,1],endpoints[k,1]]
	    d12 <- dis[endpoints[j,1],endpoints[k,2]]
	    d21 <- dis[endpoints[j,2],endpoints[k,1]]
	    d22 <- dis[endpoints[j,2],endpoints[k,2]]
	    dmin <- min(d11,d12,d21,d22)
	    if (dmin == d21) {
		endpoints[i,1] <- endpoints[j,1]
		endpoints[i,2] <- endpoints[k,2]
	    }
	    
	    else if (dmin == d11) {
		endpoints[i,1] <- endpoints[j,2]
		endpoints[i,2] <- endpoints[k,2]
		dir[i,1] <- -1
	    }
	    else if (dmin == d12) {
		endpoints[i,1] <- endpoints[j,2]
		endpoints[i,2] <- endpoints[k,1]
		dir[i,1] <- -1
		dir[i,2] <- -1
	    }
	    else  {
		endpoints[i,1] <- endpoints[j,1]
		endpoints[i,2] <- endpoints[k,1]
		dir[i,2] <- -1}}
    }
    for (i in n:2) {
	if (dir[i,1] == -1) {
	    m <- merges[i,1]
	    if (m > 0) {
		m1 <- merges[m,1]
		merges[m,1] <- merges[m,2]
		merges[m,2] <- m1
		if (dir[m,1] == dir[m,2]) 
		dir[m,] <- -dir[m,] 
	    }}
	if (dir[i,2] == -1) {
	    m <- merges[i,2]
	    if (m > 0) {
		m1 <- merges[m,1]
		merges[m,1] <- merges[m,2]
		merges[m,2] <- m1
		if (dir[m,1] == dir[m,2]) 
		dir[m,] <- -dir[m,] 
	    }}	      
	
    }
    
    clusters <- as.list(1:n)
    for (i in 1:n) {
	j <- merges[i,1]
	k <- merges[i,2]
	if ((j < 0) && (k < 0)) 
	clusters[[i]] <- c(-j,-k)
	else if (j < 0)
	clusters[[i]] <- c(-j,clusters[[k]])
	else if (k < 0)
	clusters[[i]] <- c(clusters[[j]],-k)
	else clusters[[i]] <- c(clusters[[j]], clusters[[k]])}
    
    x1 <- x
    x1$merge <- merges
    x1$order <- clusters[[n]]
    x1
    
    
}


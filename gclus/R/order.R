# Given a list whose ith element contains the indices
# of objects in the ith cluster, returns a vector whose ith 
# element gives the cluster number of the ith object.

clus2memship <- 
function(clusters) {
    ans <- 1:length(unlist(clusters))
    i <- 1
    for (cl in clusters) {
	ans[cl] <- i
	i <- i+1
    }
    ans    
}



# Given a vector whose ith elements gives the cluster number of the
# ith object, returns a list whose ith element contains the indices
# of objects in the ith cluster
memship2clus <-
function(memship) {
    m <- sort(unique(memship))
    index <- seq(along=memship)
    sapply(m, function(g) index[memship==g],simplify=FALSE)
}
    

# This function accepts a "dist" or matrix of  scores and
# returns an approximate Robinson ordering, used for scatterplot matrices.
# 
order.single <-
function(merit,clusters=NULL) {
    if (is.null(clusters))
    order.hclust(merit, TRUE,method = "single")
    else {
	dis <- - merit
	if (is.matrix(dis)) {
	    dism <- dis
	    dis <- as.dist(dis) } 
	else 
	dism <- as.matrix(dis)
	n <- nrow(dism)
	
	if (n <= 2)
	clus <- 1:n
	else {
	    cind <- col(matrix(0,n,n))
	    cind <- cind[lower.tri(cind)]
	    rind <- row(matrix(0,n,n))
	    rind <- rind[lower.tri(rind)]
	    d <- cbind(as.vector(dis),rind,cind)
	    d <- d[sort.list(d[,1],),]
	    
	    if (is.null(clusters)) {
		memship <- 1:n
		clusters <- as.list(1:n)}
		else memship <- clus2memship(clusters)		
		
		m <- length(dis)
		for (i in 1:m) {
		    j <- memship[d[i,2]]
		    k <- memship[d[i,3]]
		    if (j!= k) {
			if (j > k) {
			    r <- j
			    j <- k
			    k <- r}
			memship[memship==k] <- j
			clusj <- clusters[[j]]
			clusk <- clusters[[k]]
			dll <- dism[clusj[1], clusk[1]]
			dlr <- dism[clusj[1], clusk[length(clusk)]] 
			drl <- dism[clusj[length(clusj)], clusk[1]] 
			drr <- dism[clusj[length(clusj)], clusk[length(clusk)]] 	
			mind <- min(dll,dlr,drl,drr)
			if (drl==mind)
			NULL
			else if (dlr==mind) {
			    clusj <-rev(clusj)
			    clusk <- rev(clusk)}
			else if (dll ==mind)
			clusj <- rev(clusj)
			else clusk <- rev(clusk)
			clusters[[j]] <- c(clusj,clusk)
		    }
		    if (length(clusters[[1]]) == n) break
		}
		clus <- clusters[[1]]}
	     clus}}
    

    
	

# This function accepts a "dist" or matrix of scores and
# returns an improved ordering, for parallel coordinate displays.

order.endlink <-
function(merit,clusters=NULL) {
    dis <- - merit
    if (is.matrix(dis)) {
	dism <- dis
	dis <- as.dist(dis) } 
    else {
	dism <- as.matrix(dis)}
    n <- nrow(dism)
    if (n <= 2)
    clus <- 1:n
    else {
	cind <- col(matrix(0,n,n))
	cind <- cind[lower.tri(cind)]
	rind <- row(matrix(0,n,n))
	rind <- rind[lower.tri(rind)]
	d <- cbind(as.vector(dis),rind,cind)
	d <- d[sort.list(d[,1],),]
	if (is.null(clusters)) {
	    memship <- 1:n
	    clusters <- as.list(1:n)}
	else memship <- clus2memship(clusters)
	m <- n*(n-1)/2
	for (i in 1:m) {
	    j <- memship[d[i,2]]
	    k <- memship[d[i,3]]
	    if (!(j == k || j == -1 || k == -1)) {
		if (j > k) {
		    r <- j
		    j <- k
		    k <- r
		}
		clusj <- clusters[[j]]
		clusk <- clusters[[k]]
		dll <- dism[clusj[1], clusk[1]]
		dlr <- dism[clusj[1], clusk[length(clusk)]] 
		drl <- dism[clusj[length(clusj)], clusk[1]] 
		drr <- dism[clusj[length(clusj)], clusk[length(clusk)]] 	
		
		mind <- min(dll,dlr,drl,drr)
		if (drl==mind)
		NULL
		else if (dlr==mind) {
		    clusj <-rev(clusj)
		    clusk <- rev(clusk)}
		else if (dll ==mind)
		clusj <- rev(clusj)
		else clusk <- rev(clusk)
		clusters[[j]] <- c(clusj,clusk)
		if (! (length(clusj) == 1))
		memship[clusj[length(clusj)]] <- -1
		if (! (length(clusk) == 1))
		memship[clusk[1]] <- -1
		memship[clusk[length(clusk)]] <- j
	    }
	    if (length(clusters[[1]]) == n) break
	}
	clus<- clusters[[1]]
    }
    
    clus
}



# This function takes a merit measure and clusters, either a vector
# giving the cluster number of the ith items, or a list whose ith element
# gives the indices of the elements in the ith cluster.
# Objects within a cluster are ordered with within.order
# and clusters are ordered with between.order.
# 
order.clusters <- function(merit,clusters,within.order = order.single, 
    between.order= order.single,...) {
    if (!is.list(clusters))
    clusters <- memship2clus(clusters)
    if (!is.matrix(merit)) 
    merit <- as.matrix(merit) 
    if (!is.null(within.order)) {
	clusl <- lapply(clusters, function(g)
	    within.order(merit[g,g],...))
	newclusl <- lapply(1:length(clusters),function(i) clusters[[i]][clusl[[i]]])
    }
    else newclusl <- clusters
    if (!is.null(between.order))
    between.order(merit,newclusl)
    else unlist(newclusl)
    
}





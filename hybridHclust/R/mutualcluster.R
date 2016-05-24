setobs <- function(d){
  # Recursively add list of observations contained in each node of d.
  # d is an object of class "dendrogram"
  # assumes labels are numbers or character strings representing numbers.
  if (attr(d,'members')>1) {
    d[[1]] <- setobs(d[[1]])
    d[[2]] <- setobs(d[[2]])
    attr(d,'obs') <- c(attr(d[[1]],'obs'),attr(d[[2]],'obs'))
  } 
  else 
    attr(d,'obs') <- as.numeric(attr(d,'label'))
  return(d)
}

setmc <- function(d,dmat){
  # Recursively determine whether each node of the dendrogram d
  # represents a mutual cluster.  
  # dmat is a distance matrix, not a distance object.
  # This is the slow part, perhaps due to passing dmat.
  if (class(dmat)=='dist') dmat <- as.matrix(dmat)
  if (attr(d,'members')>1) {
    d[[1]] <- setmc(d[[1]],dmat)
    d[[2]] <- setmc(d[[2]],dmat)
    obs <- attr(d,'obs')
    # First condition below will not check the root node; although it
    # is a mc we ignore it.
    if ((length(obs)<nrow(dmat))  && (max(dmat[obs,obs]) < min(dmat[obs,-obs])))
        attr(d,'edgetext') <- ' '
  }
  return(d)
}

mutualCluster <- function(x=NULL,distances=NULL,method='average',plot=FALSE){
	if ((is.null(x) & is.null(distances)) | (!is.null(x) &
!is.null(distances)))
		stop('Exactly one of x and distances must be provided')
	if (is.null(distances)) distances <- dist(x)
	if (class(distances)!='dist') distances <- as.dist(distances)
	if (!is.null(attr(distances,'Labels'))) attr(distances,'Labels') <-
NULL
	h <- hclust(distances,method)
	themc <- getmc(h,distances,plot)
	class(themc) <- 'mutualCluster'
	dvec <- vector('list',length(themc))
	dmat <- as.matrix(distances)
	for(i in 1:length(themc)){
		if (length(themc[[i]])==2)
			dvec[[i]] <- dmat[themc[[i]][1],themc[[i]][2],drop=F]
		else {
			dvec[[i]] <- dmat[themc[[i]],themc[[i]]]
			dvec[[i]] <- as.dist(dvec[[i]])
		}
	}
	attr(themc,'distances') <- dvec
	return(themc)
}

print.mutualCluster <- function(x,...){
	for (i in 1:length(x)) cat(i,':',x[[i]],'\n')
}

get.distances <- function(x){
	print(attr(x,'distances'))
}

getmc <- function(h,dmat,plot=FALSE){
# calculate the mutual clusters using a bottom-up clustering
# optionally displaying the results.
	if (!(h$method %in% c('single','complete','average'))) stop('hclust
object must be one of single/complete/average')
	d1 <- setmc(setobs(as.dendrogram(h)),dmat)
	if (plot) plot(d1)
	return(getmc.from.dendrogram(d1))
}

getmc.from.dendrogram <- function(d){
  # Extract from a dendrogram a list of mutual clusters.
  # The dendrogram must have had "getmc" run on it already.
  if (!is.null(attr(d,'edgetext'))) 
    return(list(sort(attr(d,'obs'))))
  else {
    if (attr(d,'members')>1)
      return(c(getmc.from.dendrogram(d[[1]]),getmc.from.dendrogram(d[[2]])))
    else 
      return(NULL)
  }
}


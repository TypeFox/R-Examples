distmc <- function(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
{
    if (!is.na(pmatch(method, "euclidian"))) 
        method <- "euclidean"
    METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
        "binary", "minkowski")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
    if(!is.numeric(as.matrix(x)))
	  	stop('the data contains non-numeric data type')
    N <- nrow(x <- as.matrix(x))
    d <- .C("Rdistance", x = as.double(x), nr = N, nc = ncol(x), 
        d = double(N * (N - 1)/2), diag = as.integer(FALSE), 
        method = as.integer(method), p = as.double(p), 
        #DUP = FALSE, 
        NAOK = TRUE, PACKAGE = "Rlof")$d
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1L]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    if (method == 6) 
        attr(d, "p") <- p
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}

lof <- function(data, k, cores = NULL, ...)
{
	
  if(is.null(k))
  	stop('k is missing')
  
  if(!is.numeric(k))
  	stop('k is not numeric')
  	
  if(!is.numeric(cores) && !is.null(cores))
  	stop('cores is not numeric')

  data <- as.matrix(data)
  
  if(!is.numeric(data))
  	stop('the data contains non-numeric data type')
  
  v.k<-as.integer(k)
  
  if(max(v.k) >= dim(data)[1])
  	stop('the maximum k value has to be less than the length of the data')
  

# obtain the k nearest neighbors and their distance from each observation
  distdata <- f.dist.to.knn(data,max(v.k), cores, ...)
  
  p <- dim(distdata)[2L]
 
  # calculate the local reachability density for each observation in data

  dist.start <- as.integer((dim(distdata)[1])/2)
  dist.end <- dim(distdata)[1]
  ik <- numeric()

  registerDoParallel(cores=cores)
  m.lof <- foreach(ik = v.k, .combine=cbind) %dopar% 
  {
  	lrddata <- f.reachability(distdata,ik)
	v.lof <- rep(0,p)

  # compute the local outlier factor of each observation in data
  	for (i in 1:p)
  	{	
    	nneigh <- sum(!is.na(distdata[c((dist.start+1):dist.end),i]) & (distdata[c((dist.start+1):dist.end),i] <= distdata[(dist.start + ik),i]))
    	v.lof[i] <- sum(lrddata[distdata[(1:nneigh),i]]/lrddata[i])/nneigh
  	}
	v.lof
  # return lof, a vector with the local outlier factor of each observation
  }
  if (length(v.k) >1)
	colnames(m.lof) <- v.k
  return(m.lof)
}

f.dist.to.knn <- function(dataset,neighbors,cores,...)
{	
	m.dist <- as.matrix(distmc(dataset, ...))
	num.col <- dim(m.dist)[2]
	
	l.knndist<-lapply(c(1:num.col),function(i)
	{
		order.x <- order(m.dist[,i])
		kdist<-m.dist[,i][order.x[neighbors+1]]
		numnei <- sum(m.dist[,i] <= kdist)
		data.frame(v.order = order.x[2:numnei], v.dist = m.dist[,i][order.x[2:numnei]])
	})
	
	rm(m.dist)
	maxnum <- max(unlist(lapply(l.knndist,function(x){dim(x)[1]})))
	
	registerDoParallel(cores=cores)
	i <- numeric()
	knndist <- foreach(i = 1:num.col, .combine=cbind) %dopar%
	{
		len <- dim(l.knndist[[i]])[1]
		c(l.knndist[[i]]$v.order,rep(NA,(maxnum-len)),l.knndist[[i]]$v.dist,rep(NA,(maxnum-len)))
	}
	knndist
}

f.reachability <- function(distdata,k)
{
  p <- dim(distdata)[2]
  lrd <- rep(0,p)

  dist.start <- as.integer((dim(distdata)[1])/2)
  dist.end <- dim(distdata)[1]
    
  for (i in 1:p)
  {
    # compare the k-distance from each observation to its kth neighbor
    # to the actual distance between each observation and its neighbors
    numneigh <- sum(!is.na(distdata[c((dist.start+1):dist.end),i]) & (distdata[c((dist.start+1):dist.end),i] <= distdata[(dist.start + k),i]))
    j <- c(1:numneigh)
    temp <- rbind(distdata[dist.start+k,distdata[j,i]],distdata[dist.start+j,i])
    #calculate reachability
    reach <- 1/(sum(apply(temp,2,max))/numneigh)
    lrd[i] <- reach
  }
  lrd
}
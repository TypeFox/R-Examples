ccd.nonsequential <- function(data,m=1,method=NULL)
{
   r <- rep(0,nrow(data))
   stats <- rep(0,nrow(data))
	walks <- matrix(0,ncol=nrow(data),nrow=nrow(data))
	fs <- matrix(0,ncol=nrow(data),nrow=nrow(data))
	rx <- matrix(0,ncol=nrow(data),nrow=nrow(data))
	for(i in 1:length(r)){
	   y <- data[i,,drop=FALSE]
	   d <- as.vector(proxy::dist(data,y,method=method))
	   od <- sort(d)
	   f <- (1:nrow(data))/nrow(data)
	   dif <- f-m*(od/max(od))^2
	   r[i] <- od[which.max(dif)]
	   stats[i] <- max(dif)
	   rx[i,] <- od
	   walks[i,] <- f
	   fs[i,] <- (od/max(od))^2
	}
   n <- nrow(data)
   A <- matrix(0,nrow=n,ncol=n)
   for(i in 1:n){
      A[i,] <- proxy::dist(data[i,,drop=FALSE],data,method=method)<r[i] 
   }
	diag(A) <- 0
	out <- graph.adjacency(A,mode="Directed")
	out$R <- r
	out$stats <- stats
	out$layout <- data
	out$walks <- walks
	out$fs <- fs
	out$m <- m
	out$method <- method
	out
}

ccd.sequential <- function(data,m=1,alpha=0.05,method=NULL)
{
   n <- nrow(data)
   r <- rep(0,n)
   ks <- sqrt(qchisq(1-alpha,2)/(4*n))
   stats <- rep(0,n)
	walks <- matrix(0,ncol=n,nrow=n)
	fs <- matrix(0,ncol=n,nrow=n)
	rx <- matrix(0,ncol=n,nrow=n)
	for(i in 1:n){
	   y <- data[i,,drop=FALSE]
	   d <- as.vector(proxy::dist(data,y,method=method))
	   od <- sort(d)
	   f <- (1:nrow(data))/nrow(data)
	   dif <- f-m*(od/max(od))^2
      # find first one less than -ks
      a <- match(TRUE,c(0,dif)[1:n]< -ks,nomatch=n)
      r[i] <- od[which.max(dif[1:a])]
	   stats[i] <- max(dif[1:a])
	   rx[i,] <- od
	   walks[i,] <- f
	   fs[i,] <- (od/max(od))^2
	}
   A <- matrix(0,nrow=n,ncol=n)
   for(i in 1:n){
      A[i,] <- proxy::dist(data[i,,drop=FALSE],data,method=method)<r[i] 
   }
	diag(A) <- 0
	out <- graph.adjacency(A,mode="Directed")
	out$R <- r
	out$stats <- stats
	out$layout <- data
	out$walks <- walks
	out$fs <- fs
	out$m <- m
	out$alpha <- alpha
	out$method <- method
	out
}

ccd <- function(data,m=1,alpha=0.05,sequential=TRUE,method=NULL)
{
   if(sequential)
	   z <- ccd.sequential(data,m,alpha,method=method)
   else
	   z <- ccd.nonsequential(data,m,method=method)
	class(z) <- c("ccd",class(z))
   z
}

plot.ccd <- function(x,...)
{
	g <- x
	class(g) <- c("cccd","igraph")
	plot(g,...)
}

rng <- function(x=NULL,dx=NULL,r=1, method=NULL,usedeldir=TRUE,open=TRUE,k=NA,
                algorithm='cover_tree')
{
	if(is.na(k)){
		if(is.null(dx)) {
		  if(is.null(x)) stop("One of x or dx must be given.")
		  dx <- as.matrix(proxy::dist(x,method=method))
		} else {
			usedeldir <- FALSE
		}
		n <- nrow(dx)
		A <- matrix(0,nrow=n,ncol=n)
		if(is.vector(x)) x <- matrix(x,ncol=1)
		if(usedeldir && ncol(x)==2){
		  del <- deldir(x[,1],x[,2])
		  for(edge in 1:nrow(del$delsgs)){
			  i <- del$delsgs[edge,5]
			  j <- del$delsgs[edge,6]
			 d <- min(apply(cbind(dx[i,-c(i,j)],dx[j,-c(i,j)]),1,max))
			 if(open){
				 if(r*dx[i,j] < d){
					 A[i,j] <- 1
					 A[j,i] <- 1
				 }
			 } else {
				 if(r*dx[i,j] <= d){
					 A[i,j] <- 1
					 A[j,i] <- 1
				 }
			 }
		  }
		}
		else{
		  diag(dx) <- Inf
		  for(i in 1:n){
			 for(j in setdiff(1:n,i)){
				 d <- min(apply(cbind(dx[i,-c(i,j)],dx[j,-c(i,j)]),1,max))
				 if(open){
					 if(r*dx[i,j] < d){
						 A[i,j] <- 1
						 A[j,i] <- 1
					 }
				 } else {
					 if(r*dx[i,j] <= d){
						 A[i,j] <- 1
						 A[j,i] <- 1
					 }
				 }
			 }
		  }
		}
		diag(A) <- 0
		out <- graph.adjacency(A,mode="undirected")
	} else {
	  if(is.null(x)) stop("x must not be null")
	  n <- nrow(x)
	  k <- min(k,n-1)
	  dx <- get.knn(x,k=k,algorithm=algorithm)
	  edges <- NULL
	  for(i in 1:n){
		  i.indices <- dx$nn.index[i,]
		  i.dists <- dx$nn.dist[i,]
		  for(j in 1:k){
			  rd <- r*i.dists[j]/2
			  j.indices <- dx$nn.index[i.indices[j],]
			  j.dists <- dx$nn.dist[i.indices[j],]
			  rd <- r*i.dists[j]
			  S <- setdiff(intersect(i.indices,j.indices),c(i,i.indices[j]))
			  if(length(S)>0){
				  d <- Inf
				  for(si in S){
					  a <- which(i.indices == si)
					  b <- which(j.indices == si)
					  d <- min(d,max(i.dists[a],j.dists[b]))
				  }
				  if(rd < d){
					  edges <- cbind(edges,c(i,i.indices[j]))
				  }
			  }
		  }
	  }
	  out <- simplify(graph(edges,n=n,directed=FALSE))
	}
	if(!is.null(x)){
		out$layout <- x
	}
	out$r <- r
   out
}


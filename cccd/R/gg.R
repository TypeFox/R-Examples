gg <- function(x,r=1,method=NULL,usedeldir=TRUE,open=TRUE,k=NA,algorithm='cover_tree')
{
	n <- nrow(x)
	if(is.na(k)){
		dx <- as.matrix(proxy::dist(x,method=method))
		A <- matrix(0,nrow=n,ncol=n)
		if(is.vector(x)) x <- matrix(x,ncol=1)
		if(usedeldir && (method=="Euclidean" || is.null(method)) && (ncol(x)==2)){
		  del <- deldir(x[,1],x[,2])
		  for(edge in 1:nrow(del$delsgs)){
			 i <- del$delsgs[edge,5]
			 j <- del$delsgs[edge,6]
			 d1 <- r*dx[i,j]/2
			 d <- proxy::dist((x[i,,drop=FALSE]+x[j,,drop=FALSE])/2,x,
			          method=method)
			 d[i] <- Inf
			 d[j] <- Inf
			 if(open){
				 if(!any(d<d1)){
					 A[i,j] <- 1
					 A[j,i] <- 1
				 }
			 } else {
				 if(!any(d<=d1)){
					 A[i,j] <- 1
					 A[j,i] <- 1
				 }
			 }
		  }
		} else {
			for(i in 1:n){
			  for(j in setdiff(1:n,i)){
				 d1 <- r*dx[i,j]/2
				 d1 <- round(d1,10)
				 d <- proxy::dist((x[i,,drop=FALSE]+x[j,,drop=FALSE])/2,x,
				          method=method)
				 d <- round(d,10)
				 d[i] <- Inf
				 d[j] <- Inf
				 if(open){
					 if(!any(d<d1)){
						 A[i,j] <- 1
						 A[j,i] <- 1
					 }
				 } else {
					 if(!any(d<=d1)){
						 A[i,j] <- 1
						 A[j,i] <- 1
					 }
				 }
			 }
		  } 
   } 
	diag(A) <- 0
	g <- graph.adjacency(A,mode="undirected")
  } else {
	  k <- min(k,n-1)
	  dx <- get.knn(x,k=k,algorithm=algorithm)
	  edges <- NULL
	  for(i in 1:n){
		  i.indices <- dx$nn.index[i,]
		  i.dists <- dx$nn.dist[i,]
		  for(j in 1:k){
			  rd <- r*i.dists[j]/2
			  j.indices <- dx$nn.index[i.indices[j],]
			  S <- setdiff(intersect(i.indices,j.indices),c(i,i.indices[j]))
			  if(length(S)>0){
				  X <- (x[i,,drop=FALSE]+x[i.indices[j],,drop=FALSE])/2
				  Y <- x[S,,drop=FALSE]
				  d <- min(as.vector(proxy::dist(X,Y)))
				  if(d>rd){
					  edges <- cbind(edges,c(i,i.indices[j]))
				  }
			  }
		  }
	  }
	  g <- simplify(graph(edges,n=n,directed=FALSE))
  }
	g$layout <- x
	g$r <- r
   g
}


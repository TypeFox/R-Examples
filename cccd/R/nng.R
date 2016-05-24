
nng <- function(x=NULL,dx=NULL,k=1,mutual=FALSE,method=NULL,
                use.fnn=FALSE,algorithm='cover_tree')
{
	if(use.fnn){
	   if(is.null(x)) stop("x must not be null. try use.fnn=FALSE")
		dx <- get.knn(x,k=k,algorithm=algorithm)
		edges <- matrix(unlist(sapply(1:nrow(x),function(i) {
					 rbind(rep(i,k),dx$nn.index[i,])
			 })),nrow=2)
		n <- nrow(x)
		if(mutual){
			a <- apply(edges,2,sort)
			b <- which(duplicated(a,MARGIN=2))
			if(length(b)==0){
			   out <- graph.empty(n,directed=FALSE)
			} else {
				out <- graph(edges[,b,drop=FALSE],n=n,directed=FALSE)
			}
		} else {
			out <- graph(edges,n=n,directed=TRUE)
		}
	} else {
		if(is.null(dx)) {
		  if(is.null(x)) stop("one of x or dx must be given")
		  dx <- as.matrix(proxy::dist(x,method=method))
		}
		n <- nrow(dx)
		A <- matrix(0,nrow=n,ncol=n)
		for(i in 1:n){
		  d <- sort(dx[i,])[-1]
		  A[i,dx[i,]<=d[k]] <- 1
		}
		diag(A) <- 0
		if(mutual){
			for(i in 1:n){
			  A[i,] <- A[i,] & A[,i]
			  A[,i] <- A[i,]
		  }
		}
		if(mutual)
			out <- graph.adjacency(A,mode="undirected")
		else
			out <- graph.adjacency(A,mode="directed")
	}
	out$k <- k
	out$mutual <- mutual
	if(!is.null(x)){
		out$layout <- x
	}
   out
}


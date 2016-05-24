vote <- function(x)
{
   cls <- unique(x)
	best <- 0
	k <- 0
	for(class in cls){
	   a <- sum(x==class)
	   if(a>best){
		   best <- sum(cls==class)
			k <- class
		}
		else if(a==best){
		   k <- c(k,class)
		}
	}
   k
}

prune <- function(x,classes,prox="Gabriel",ignore.ties=TRUE,...)
{
	MODES <- c("Gabriel","Relative Neighborhood","k-Nearest Neighbor","Minimum Spanning Tree")
	tmp <- charmatch(prox,MODES)
	if(is.null(tmp)){
		 stop("invalid proximity graph or proximity graph not recognized")
	}
	else if(is.na(tmp)){
		 stop("invalid proximity graph or proximity graph not recognized")
	}
	else if(tmp==0){
		 stop("ambiguous proximity graph: retry with more characters")
	}
	else if(tmp<1 || tmp>length(MODES)){
		 stop("invalid proximity graph or proximity graph not recognized")
	}
	mode <- MODES[tmp]
	if(mode=="Gabriel")
		g <- gg(x,...)
	else if(mode=="Relative Neighborhood")
		g <- rng(x,...)
	else if(mode=="k-Nearest Neighbor")
		g <- nng(x,...)
	else if(mode=="Minimum Spanning Tree"){
		D <- as.matrix(proxy::dist(x))
		n <- vcount(g)
		A <- matrix(1,nrow=n,ncol=n)
		diag(A) <- 0
		h <- graph.adjacency(A,mode="undirected")
		w <- rep(0,choose(n,2))
		k <- 1
		for(i in 1:(n-1)){
			for(j in (i+1):n){
					w[k] <- D[i,j]
					k <- k+1
			}
		}
		g <- minimum.spanning.tree(h,weights=w,...)
	}
	n <- vcount(g)
	v <- NULL
	for(i in 1:n){
	   a <- setdiff(neighborhood(g,order=1,nodes=i),i)
		w <- vote(classes[a]) 
		if(ignore.ties){
			if(any(classes[i] %in% w)) v <- c(v,i)
		}
		else {
			if(length(w)==1)
				if(classes[i] == w) v <- c(v,i)
		}
	}
	g$layout <- x
	list(x=x[v,],v=v,graph=g)
}

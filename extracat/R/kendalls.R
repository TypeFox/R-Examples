kendalls <- function(x){
	cs <- colSums(x)
	rs <- rowSums(x)
	css <- rev(cumsum(rev(cs)))
	rss <- rev(cumsum(rev(rs)))
	n <- nrow(x)
	m <- ncol(x)
	N <- sum(x)
	xx <- sum(cs*(css-cs)[1:m])
	y <- sum(rs*(rss-rs)[1:n])
	
	
	storage.mode(x) <- "integer"
	x2 <- x[,m:1]
	
	dims <- as.integer(c(n,m))
	crt1 <- .Call("classcrit",x,dims,as.integer(0))
	crt2 <- .Call("classcrit",x2,dims,as.integer(0))
	
	tau <- (crt1-crt2)/sqrt(max(xx,1))/sqrt(max(1,y))
#scrt <- crt1/x/y*N^2
#return(c(tau,crt1,scrt,n,m,N))
	return(tau)
}


kendallsMV <- function(x){
	
	nd <- length(dim(x))
	ms <- list()
	mss <- list()
	for(i in 1:nd){
		ms[[i]] <- apply(x,i,sum)
		mss[[i]] <- rev(cumsum(rev(ms[[i]])))
	}
	N <- sum(x)
		storage.mode(x) <- "integer"
	crt1 <- .Call("classcrit",x,as.integer(dim(x)),as.integer(0))
	crt2 <- allpairs(x)-crt1
	
	xx <- mapply(function(y,z){
		 max(1,sum(y*(z-y)[1:length(y)]))^(1/nd)
	},y = ms,z=mss)
	
	tau <- (crt1-crt2)/prod(xx)
	return(tau)
}
juggle1 <- function(data,classes,sampled=TRUE,
                   num=100,sample.proportion=0.1,method=NULL)
{
	if(sample.proportion<1.0) sampled<-TRUE
	replace <- sample.proportion>=1
	cls <- 1:length(unique(classes))
	if(is.null(data)) stop("data must be given")
	if(is.matrix(data)){
		d <- ncol(data)
	}
	else d <- 1
	dxx <- as.matrix(proxy::dist(data,method=method))
	n <- nrow(dxx)
	S <- list(0)
	R <- list(0)
	for(i in 1:num){
		S[[i]] <- list(0)
		R[[i]] <- list(0)
		xindex <- list(0)
		for(j in cls){
			xindex[[j]] <- which(classes==j)
			if(sampled){
				nx <- sample.proportion*length(xindex[[j]])
				xindex[[j]] <- sort(unique(sample(xindex[[j]],nx,replace=replace)))
			}
		}
		for(j in cls){
			DXX <- dxx[xindex[[j]],xindex[[j]],drop=FALSE]
			ind <- setdiff(unlist(xindex),xindex[[j]])
			DYX <- dxx[ind,xindex[[j]],drop=FALSE]
			rx <- apply(DYX,2,min)
			A <- matrix(as.integer(DXX<rx),nrow=nrow(DXX))
			diag(A) <- 0
			Gx <- graph.adjacency(A,mode="directed")
			sx <- dominate.random.sample(Gx)
			S[[i]][[j]] <- xindex[[j]][sx]
			R[[i]][[j]] <- rx[sx]
		}
	}
	list(S=S,R=R,dimension=d,method=method)
}

juggle1.classify <- function(data,J,tdata,indices)
{
	if(missing(indices)) indices <- 1:length(J$S)
	n <- length(indices)
	if(any(indices>length(J$S))) stop("invalid indices")
	if(J$dimension>1)
		N <- nrow(data)
	else N <- length(data)
	nc <- length(J$S[[1]])
	out <- matrix(0,nrow=N,ncol=nc)
	inds <- sort(unique(unlist(J$S)))
	d <- as.matrix(proxy::dist(data,tdata[inds,,drop=FALSE],method=J$method))
	foo <- matrix(0,nrow=nc,ncol=N)
	for(i in indices){
		for(j in 1:nc){
			if(length(J$S[[i]][[j]])>1)
				if(N==1)
					foo[j,] <- min(d[,match(J$S[[i]][[j]],inds)]/J$R[[i]][[j]])
				else
					foo[j,] <- apply(t(t(d[,match(J$S[[i]][[j]],inds)])/J$R[[i]][[j]]),1,min)
			else
				foo[j,] <- d[,match(J$S[[i]][[j]],inds)]/J$R[[i]][[j]]
		}
		a <- apply(foo,2,which.min)
		for(j in 1:nc){
			out[,j] <- out[,j]+1*(a==j)
		}
	}
	out/n
}

juggle2 <- function(data,classes,sampled=TRUE,
                   num=100,sample.proportion=0.1,
						 k=2,method=NULL)
{
	if(sample.proportion<1.0) sampled<-TRUE
	replace <- sample.proportion>=1
	cls <- 1:length(unique(classes))
	if(is.null(data)) stop("data must be given")
	if(is.matrix(data)){
		d <- ncol(data)
	}
	else d <- 1
	if(k>=d) stop("use juggle if no dimension sampling is to be performed")
	n <- nrow(data)
	S <- list(0)
	R <- list(0)
	vars <- list(0)
	for(i in 1:num){
		S[[i]] <- list(0)
		R[[i]] <- list(0)
		xindex <- list(0)
		for(j in cls){
			xindex[[j]] <- which(classes==j)
			if(sampled){
				nx <- sample.proportion*length(xindex[[j]])
				xindex[[j]] <- sort(unique(sample(xindex[[j]],nx,replace=replace)))
			}
		}
		if(k>0) K<-k
		else{
		   K <- rbinom(1,d,.5)
		}
		vars[[i]] <- sample(1:d,K)
		dxx <- as.matrix(proxy::dist(data[,vars[[i]],drop=FALSE],method=method))
		for(j in cls){
			DXX <- dxx[xindex[[j]],xindex[[j]]]
			ind <- setdiff(unlist(xindex),xindex[[j]])
			DYX <- dxx[ind,xindex[[j]]]
			rx <- apply(DYX,2,min)
			A <- matrix(as.integer(DXX<rx),nrow=nrow(DXX))
			diag(A) <- 0
			Gx <- graph.adjacency(A,mode="directed")
			sx <- dominate.random.sample(Gx)
			S[[i]][[j]] <- xindex[[j]][sx]
			R[[i]][[j]] <- rx[sx]
		}
	}
	list(S=S,R=R,dimension=d,vars=vars,method=method)
}

juggle2.classify <- function(data,J,tdata,indices)
{
	if(missing(indices)) indices <- 1:length(J$S)
	n <- length(indices)
	if(any(indices>length(J$S))) stop("invalid indices")
	if(J$dimension>1)
		N <- nrow(data)
	else N <- length(data)
	nc <- length(J$S[[1]])
	out <- matrix(0,nrow=N,ncol=nc)
	inds <- sort(unique(unlist(J$S)))
	foo <- matrix(0,nrow=nc,ncol=N)
	for(i in indices){
		d <- as.matrix(proxy::dist(data[,J$vars[[i]],drop=FALSE],tdata[inds,J$vars[[i]],drop=FALSE],method=J$method))
		for(j in 1:nc){
			if(length(J$S[[i]][[j]])>1)
				if(N==1)
					foo[j,] <- min(d[,match(J$S[[i]][[j]],inds)]/J$R[[i]][[j]])
				else
					foo[j,] <- apply(t(t(d[,match(J$S[[i]][[j]],inds)])/J$R[[i]][[j]]),1,min)
			else
				foo[j,] <- d[,match(J$S[[i]][[j]],inds)]/J$R[[i]][[j]]
		}
		a <- apply(foo,2,which.min)
		for(j in 1:nc){
			out[,j] <- out[,j]+1*(a==j)
		}
	}
	out/n
}


juggle <- function(data,classes,sampled=TRUE,sample.dim=FALSE,
                   num=100,sample.proportion=0.1,
						 k=2,method=NULL)
{
   if(sample.dim) a <- juggle2(data,classes,sampled,num,
	                            sample.proportion,k,method=method)
   else a <- juggle1(data,classes,sampled,num,sample.proportion,method=method) 
	a
}

juggle.classify <- function(data,J,tdata,indices)
{
   if(is.null(J$vars)) a <- juggle1.classify(data,J,tdata,indices)
   else a <- juggle2.classify(data,J,tdata,indices)
	a
}

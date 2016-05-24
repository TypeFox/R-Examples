cccd <- function(x=NULL,y=NULL,dxx=NULL,dyx=NULL,method=NULL,k=NA,algorithm='cover_tree')
{
	if(is.null(dyx) && (is.null(method) || method=='euclidean')){
	   dyx <- get.knnx(y,x,k=1,algorithm=algorithm)
		R <- dyx$nn.dist[,1]
	} else if(is.null(dyx)){
		dyx <- as.matrix(proxy::dist(y,x,method=method))
		R <- apply(dyx,2,min)
	} else {
		R <- apply(dyx,2,min)
	}
	if(is.na(k)){
		if(is.null(dxx) | is.null(dyx)){
			if(is.null(x) | is.null(y)) stop("either x,y or dxx,dyx must be given")
			dxx <- as.matrix(proxy::dist(x,method=method))
		}
		M <- matrix(as.integer(dxx<R),length(R))
		diag(M) <- 0
		g <- graph.adjacency(M,mode="directed")
	} else {
		if(is.null(x) || is.null(y)) stop("x and y must not be null")
		k <- min(k,nrow(y))
	   dyx <- get.knnx(y,x,k=1,algorithm=algorithm)
	   dxx <- get.knn(x,k=k,algorithm=algorithm)
		R <- dyx$nn.dist[,1]
		out <- unlist(sapply(1:nrow(x), function(i) {
						a <- which(dxx$nn.dist[i,]<R[i])
						if(length(a)==0) return(NULL)
		            rbind(rep(i, length(a)), dxx$nn.index[i,a])
					}))
		if(is.null(out)){
			g <- graph.empty(n=nrow(x),directed=TRUE)
		} else {
			edges <- matrix(out,nrow=2)
			g <- graph(edges,n=nrow(x),directed=TRUE)
		}
	}
	g$R <- R
	if(is.vector(x) || (dim(x)==1)) {
		x <- cbind(x,rep(0,length(x)))
		y <- cbind(y,rep(0,length(y)))
	}
	g$layout <- x
	g$Y <- y
	g$method <- method
	class(g) <- c("cccd",class(g))
	g
}

cccd.rw <- function(x=NULL,y=NULL,dxx=NULL,dyx=NULL,method=NULL,m=1,d=2)
{
   if(is.null(dxx) | is.null(dyx)){
		if(is.null(x) | is.null(y)) stop("either x,y or dxx,dyx must be given")
      dyx <- as.matrix(proxy::dist(y,x,method=method))
      dxx <- as.matrix(proxy::dist(x,method=method))
		d <- ncol(x)
	} 
   R <- rep(0,nrow(dxx))
	nx <- nrow(dxx)
	ny <- nrow(dyx)
	for(i in 1:nx){
		o <- order(c(dxx[i,],dyx[,i]))
	   rw <- cumsum(c(rep(1/nx,nx),rep(-1/ny,ny))[o])
		r <- sort(dxx[i,])
		R[i] <- r[which.max(rw[1:nx]-m*r^d)]
	}
	M <- matrix(as.integer(dxx<R),length(R))
	diag(M) <- 0
	g <- graph.adjacency(M,mode="directed")
	g$R <- R
	g$layout <- x
	g$Y <- y
	g$method <- method
	class(g) <- c("cccd",class(g))
	g
}

plot.cccd <- function(x,...,
                     plot.circles=FALSE,dominate.only=FALSE,D=NULL,
							vertex.size=2,vertex.label=NA,
							vertex.color="SkyBlue2",dom.color="Blue",
                     ypch=20,
							ycex=1.5,ycol=2,
							use.circle.radii=FALSE,
							balls=FALSE,
							ball.color=gray(.8),
							square=FALSE,
							xlim,ylim)

{
	g <- x
	class(g) <- "igraph"
	if(balls) plot.circles <- TRUE
	x <- g$layout
	n <- nrow(x)
	y <- g$Y
	if(is.null(y)){
		if(missing(xlim)){
			xlim <- range(x[,1])
		}
		if(missing(ylim)){
			ylim <- range(x[,2])
		}
	}
	else {
		if(missing(xlim)){
			xlim <- range(c(x[,1],y[,1]))
		}
		if(missing(ylim)){
			ylim <- range(c(x[,2],y[,2]))
		}
	}
	if(is.null(D)) D <- dominate(g)
	col <- rep(vertex.color,n)
	col[D] <- dom.color
	vertex.color <- col
	r <- g$R
	if(use.circle.radii){
		 xlim <- range(c(xlim[2]+r,xlim[1]-r))
		 ylim <- range(c(ylim[2]+r,ylim[1]-r))
	}
	if(square){
	   xlim <- range(c(xlim,ylim))
	   ylim <- xlim
	}
	plot(g,xlim=xlim,ylim=ylim,vertex.size=vertex.size,rescale=FALSE,
	           vertex.label=vertex.label,
				  vertex.color=vertex.color,...)
	if(plot.circles){
		col <- rep(ifelse(balls,ball.color,vertex.color),n)
		if(dominate.only){
		   col[-D] <- NA
		}
		if(balls){
			draw.circle(x[,1],x[,2],r,border=vertex.color,col=col)
			plot(g,xlim=xlim,ylim=ylim,vertex.size=vertex.size,rescale=FALSE,
				  vertex.label=vertex.label,
				  vertex.color=vertex.color,add=TRUE,...)
		} else {
			draw.circle(x[,1],x[,2],r,border=col)
		}
	}
	if(!is.null(y)){
		points(y,pch=ypch,col=ycol,cex=ycex)
	}
}

cccd.classifier <- function(x,y,dom.method='greedy',proportion=1,...)
{
	if(missing(y)){
	   if(is.list(x)){
		   y <- x$y
			x <- x$x
			if(is.null(x) | is.null(y))
			   stop("must provide either x and y or a list with attributes x and y")
		}
		else
			stop("must provide either x and y or a list with attributes x and y")
	}
   Gx <- cccd(x,y,...)
	Gy <- cccd(y,x,...)
	Dx <- dominate(Gx,method=dom.method,proportion=proportion)
	Dy <- dominate(Gy,method=dom.method,proportion=proportion)
	h <- list(Rx=Gx$R[Dx],Ry=Gy$R[Dy],Cx=matrix(x[Dx,],ncol=ncol(x)),Cy=matrix(y[Dy,],ncol=ncol(y)))
	class(h) <- "cccdClassifier"
	h
}

plot.cccdClassifier <- function(x,...,xcol=1,ycol=2,xpch=20,ypch=xpch,
                                balls=FALSE,add=FALSE)
{
   data <- rbind(x$Cx,x$Cy)
	R <- c(x$Rx,x$Ry)
	cols <- c(rep(xcol,length(x$Rx)),rep(ycol,length(x$Ry)))
	pchs <- c(rep(xpch,length(x$Rx)),rep(ypch,length(x$Ry)))
	if(add){
		points(data,col=cols,pch=pchs,...)
	} else {
		plot(data,col=cols,pch=pchs,...)
	}
	if(is.character(balls) || balls){
		if(is.character(balls)){
		   if('x' %in% balls){
				draw.circle(x$Cx[,1],x$Cx[,2],x$Rx,border=cols,col=xcol)
			}
		   if('y' %in% balls){
				draw.circle(x$Cy[,1],x$Cy[,2],x$Ry,border=cols,col=ycol)
			}
			draw.circle(data[,1],data[,2],R,border=cols)
		} else {
			draw.circle(data[,1],data[,2],R,border=cols,col=cols)
		}
	} else {
		draw.circle(data[,1],data[,2],R,border=cols)
	}
}

cccd.classifier.rw <- function(x,y,m=1,d=2)
{
	if(missing(y)){
	   if(is.list(x)){
		   y <- x$y
			x <- x$x
			if(is.null(x) | is.null(y))
			   stop("must provide either x and y or a list with attributes x and y")
		}
		else
			stop("must provide either x and y or a list with attributes x and y")
	}
   Gx <- cccd.rw(x,y,m=m,d=d)
	Gy <- cccd.rw(y,x,m=m,d=d)
	Dx <- dominate(Gx)
	Dy <- dominate(Gy)
	h <- list(Rx=Gx$R[Dx],Ry=Gy$R[Dy],
	          Cx=matrix(x[Dx,],ncol=ncol(x)),Cy=matrix(y[Dy,],ncol=ncol(y)))
	class(h) <- "cccdClassifier"
	h
}

cccd.classify <- function(data,C,method=NULL)
{
	dx <- apply(t(t(as.matrix(proxy::dist(data,C$Cx,method=method)))/C$Rx),1,min)
	dy <- apply(t(t(as.matrix(proxy::dist(data,C$Cy,method=method)))/C$Ry),1,min)
	dx<dy
}

cccd.multiclass.classifier <- function(data,classes,
                                       dom.method='greedy',proportion=1,...)
{
	cls <- unique(classes)
	nc <- length(cls)
	G <- list(0)
	D <- list(0)
	C <- list(0)
	R <- list(0)
	for(i in 1:nc){
	   z <- classes==cls[i]
		x <- data[z,]
		y <- data[!z,]
		G[[i]] <- cccd(x,y,...)
		D[[i]] <- dominate(G[[i]],method=dom.method,proportion=proportion)
		C[[i]] <- matrix(x[D[[i]],],ncol=ncol(x))
		R[[i]] <- G[[i]]$R[D[[i]]]
	}
	list(G=G,D=D,C=C,R=R,classes=cls)
}

cccd.multiclass.classify <- function(data,C,method=NULL)
{
	nc <- length(C$R)
	if(is.vector(data)) data <- matrix(data,nrow=1)
	d <- matrix(0,nrow=nc,ncol=nrow(data))
	classes <- C$classes
	for(i in 1:nc){
		d[i,] <- apply(t(t(as.matrix(proxy::dist(data,C$C[[i]],method=method)))/C$R[[i]]),1,min)
	}
	z <- t(apply(d,2,function(x)x/sum(x)))
	list(probs=z,classes=classes[apply(z,1,which.min)])
}


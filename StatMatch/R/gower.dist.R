`gower.dist` <-
function (data.x, data.y=data.x, rngs=NULL, KR.corr=TRUE)
{

######
# a function to compute the Gower distance among a two vectors x and y
#
	gower.fcn <- function(x, y, rng=NULL, KR.corr=TRUE)
	{
		nx <- length(x)
		ny <- length(y)
		cx <- class(x) 
		cy <- class(y)
		delta <- matrix(1, nx, ny)
		if(!identical(cx, cy)) stop("the x and y object are of different type")
		if(is.logical(x)){
			dd <- abs(outer(X=x, Y=y, FUN="-"))
			delta[outer(x==FALSE, y==FALSE, FUN="&")] <- 0
			delta[outer(is.na(x), is.na(y), FUN="|")] <- 0
		}
		else if(is.character(x) || (is.factor(x) && !is.ordered(x)) ){
			if(is.factor(x) && !identical(levels(x), levels(y))) stop("x and y have different levels")
			xx <- c(matrix(as.character(x), nx, ny))
			yy <- c(matrix(as.character(y), nx, ny, byrow=TRUE))
			dd <- 1 - outer(x, y, FUN="==")
			delta[outer(is.na(x), is.na(y), FUN="|")] <- 0
		}	
		else if(is.ordered(x)){
			if(KR.corr){
				x <- as.numeric(x)
				y <- as.numeric(y)
				if(is.null(rng) || is.na(rng)) rng <- max(x,y) - 1
				zx <- (x-1)/rng
				zy <- (y-1)/rng
				dd <- abs(outer(X=zx, Y=zy, FUN="-"))/(max(zx,zy)-min(zx,zy) )
				delta[outer(is.na(zx), is.na(zy), FUN="|")] <- 0
			}
			else{
				x <- as.numeric(x)
				y <- as.numeric(y)
				if(is.null(rng) || is.na(rng)) rng <- max(x,y) - 1
				dd <- abs(outer(X=x, Y=y, FUN="-"))/rng
				delta[outer(is.na(x), is.na(y), FUN="|")] <- 0
			}
		}
		else{
  			if(is.null(rng) || is.na(rng)) rng <- max(x,y) - min(x,y)
	 		  dd <- abs(outer(X=x, Y=y, FUN="-"))/rng
		    delta[outer(is.na(x), is.na(y), FUN="|")] <- 0
		}	
		list(dist=dd, delta=delta)
	}
	
#########################################	
	
	if(is.null(dim(data.x)) && is.null(dim(data.y))) {
		out.gow <- gower.fcn(x=data.x, y=data.y, rng=rngs, KR.corr=KR.corr)
		out <- (out.gow$dist	* out.gow$delta) / out.gow$delta
	}
	else if(is.null(dim(data.x)) && !is.null(dim(data.y))){
		p <- ncol(data.y)
		if(length(data.x)!=p) stop("data.x should be of the same length of the no. of cols of data.y")
		num <- array(0, c(1,nrow(data.y),p))
		den <- array(0, c(1,nrow(data.y),p))
		for(k in 1:p){
			if(is.null(rngs)) rng.k <- NULL
			else rng.k <- rngs[k]
			out.gow <- gower.fcn(x=data.x[,k], y=data.y[,k], rng=rng.k, KR.corr=KR.corr)
			num[,,k] <- out.gow$dist*out.gow$delta
			den[,,k] <- out.gow$delta
		}
		out <- apply(num, c(1,2), sum, na.rm=TRUE)/apply(den, c(1,2), sum, na.rm=TRUE)
	}
	else{
		p <- ncol(data.y)
		if(ncol(data.x)!=p) stop("data.x and data.y must have the same no. of cols")
		num <- array(0, c(nrow(data.x), nrow(data.y), p))
		den <- array(0, c(nrow(data.x), nrow(data.y), p))
		for(k in 1:p){
			if(is.null(rngs)) rng.k <- NULL
			else rng.k <- rngs[k]
			out.gow <- gower.fcn(x=data.x[,k], y=data.y[,k], rng=rng.k, KR.corr=KR.corr)
			num[,,k] <- out.gow$dist*out.gow$delta
			den[,,k] <- out.gow$delta
		}
		out <- apply(num, c(1,2), sum, na.rm=TRUE)/apply(den, c(1,2), sum, na.rm=TRUE)
	}
	out
}


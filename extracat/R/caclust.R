caclust <- function(x, data = NULL, k = 3, d = NULL, link = c("complete","average","ward","single"), ... ){
	if(!inherits(x, "ca")){
		stopifnot(length(dim(x)) == 2)
		data <- as.data.frame(as.table(x))
		#require(ca)
		if (requireNamespace("ca", quietly = TRUE)) {
			x <- ca::ca(x)
		}else{
			stop("Please install package 'ca' to use this reordering function.")
		}
	}
	link <- match.arg(link)
	nr <- nrow(x$rowcoord)
	nc <- nrow(x$colcoord)
	
	d <- min(d, nr-1, nc-1)
	
	rows <- x$rowcoord[, 1:d] * rep( x$sv[1:d], each = nr)
	columns <- x$colcoord[, 1:d] * rep( x$sv[1:d], each = nc)
	
	DR <- dist(rows)
	DC <- dist(columns)
	angs <- apply(rows,1, function(z){
		apply(columns,1,function(y){
			angle(y,z)	
		})	
	})
	#print(angs)
	DM <- matrix(0,ncol=nr+nc,nrow=nr+nc)
	angs <- angs/mean(as.vector(angs))
	
	DM[1:nr,1:nr] <- as.matrix(DR/mean(as.vector(DR)))
	DM[(nr+1):(nr+nc), (nr+1):(nr+nc)] <- as.matrix(DC/mean(as.vector(DC)))
	DM[1:nr, (nr+1):(nr+nc)] <- t(angs)
	DM[(nr+1):(nr+nc), 1:nr] <- angs
	#print(DM)
	if(is.null(x$rownames)){
		x$rownames <- 1:nr	
	}
	if(is.null(x$colnames)){
		x$colnames <- 1:nc	
	}
	rown <- paste("ROW",x$rownames)
	coln <- paste("COL",x$colnames)

	rownames(DM) <- c( rown, coln )
	DM <- as.dist(DM)
	
	hc <- hclust(DM, method=link)
	#plot(hc)
	grp <- subtree(hc, k = k)$data
	ord.labs <- c(levels(data[,1]), levels(data[,2]))[hc$order]
	rr <- hc$labels[hc$order] %in% rown
	cc <- hc$labels[hc$order] %in% coln
	
	rowlabs <-  ord.labs[rr]
	collabs <-  ord.labs[cc]
	
	roword <- rank(hc$order[ rr ])
	colord <- rank(hc$order[ cc ])
	rowgrp <- grp[ hc$order ][ rr ]
	colgrp <- grp[ hc$order ][ cc ]
	
	data[,1] <- factor(data[,1], levels = rowlabs)	
	data[,2] <- factor(data[,2], levels = collabs)	
	ret <- xtabs(Freq~.,data=data)
	attr(ret,"rowgrp") <- rowgrp
	attr(ret,"colgrp") <- colgrp
	attr(ret,"tree") <- hc
	return(ret)
}

angle <- function(x,y){
	lx <- sqrt(sum( x^2 ))	
	ly <- sqrt(sum( y^2 ))
	s <- sum(x*y)/lx/ly 
	s <- min(s,1)
	ang <- acos( s )/2/pi*360
	
	return( ang )
}

plot.binaryTree <- function(x, value=TRUE, precision = 2, size=FALSE, white=TRUE, col.grid=gray.colors(max(order(tree2vec(x)))), ...)
{
	x.coord <- 0.5
	y.coord <- -0.5
	for(i in 1:x$max.s)
	{
		x.coord <- c(x.coord, seq(0,1, by=1/(2^i+1)))
		y.coord <- c(y.coord, rep(-i-0.5, 2^i))
	}
	x.coord <- x.coord[which(x.coord>0 & x.coord<1)]
	x0 <- rep(x.coord[-c((2^x$max+1):2^(x$max.s+1)-1)],each=2)
	x1 <- x.coord[-1]
	y0 <- rep(y.coord[-c((2^x$max+1):2^(x$max.s+1)-1)],each=2)
	y1 <- y.coord[-1]
	plot(0,0, ylim=c(-x$max.s-1,0), xlim=c(0,1), col=0, ylab="scales", xlab="", tcl=0, col.axis=0, ...)
	segments(x0,y0,x1,y1)
	if(white) bg = "white"
	else bg=col.grid[rank(tree2vec(x))]
	if(size)scale <- (4/max(tree2vec(x)))*tree2vec(x)
	else scale=3
	points(x.coord,y.coord, cex=scale, bg=bg, pch=21, col=1)
	if(value)text(x.coord,y.coord, round(tree2vec(x),precision), cex=scale/4)
}
#----------------
summary.binaryTree <- function(object, ...)
{
	cat("Binary Tree with ", object$max.s, "scales \n")
	print(object$T)
}
#-----------------


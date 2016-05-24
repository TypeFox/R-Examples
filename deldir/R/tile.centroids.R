tile.centroids <- function(xxx){
foo <- function (x,y) {
	x <- c(x,x[1])
	y <- c(y,y[1])
	n <- length(x)
	u <- x[-n]*y[-1]
	v <- x[-1]*y[-n]
	w <- x[-n]+x[-1]
	z <- y[-n]+y[-1]
	A <- 3*(sum(u)-sum(v))
	Cx <- sum(w*(u-v))/A
	Cy <- sum(z*(u-v))/A
	c(Cx,Cy)
}
ccc <- lapply(xxx,function(x){foo(x$x,x$y)})
ccc <- as.data.frame(matrix(unlist(ccc),ncol=2,byrow=TRUE))
names(ccc) <- c("x","y")
ccc
}

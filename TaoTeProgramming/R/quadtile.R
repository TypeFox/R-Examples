"quadtile" <-
function (dims=c(10,10), 
          color=grep("^gray", colors(), value=TRUE), seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	x <- apply(array(runif(prod(dims)), dims), 2, cumsum)
	y <- apply(array(runif(prod(dims)), dims), 2, cumsum)

	plot(1,1, type="n", xlim=range(x), ylim=range(y), xlab="",
		ylab="", axes=FALSE)
	for(i in 2:dims[1]) {
		for(j in 2:dims[2]) {
			sub <- cbind(c(i-1, i, i, i-1), c(j-1, j-1, j, j))
			polygon(x[sub], y[sub], col=safesample(color),
				border=NA)
		}
	}
}

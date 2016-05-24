"quadtilebalance" <-
function (dims=c(10,10), 
          color=grep("^gray", colors(), value=TRUE), seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	rawx <- array(runif(prod(dims)), dims)
	rawy <- array(runif(prod(dims)), dims)
	x <- apply(scale(rawx, center=FALSE, scale=colSums(rawx[-1,])), 2, 
		cumsum)
	y <- apply(scale(rawy, center=FALSE, scale=colSums(rawy[-1,])), 2, 
		cumsum)

	plot(1,1, type="n", xlim=range(x), ylim=range(y), xlab="",
		ylab="", axes=FALSE)
	for(j in 2:dims[2]) {
		for(i in 2:dims[1]) {
			sub <- cbind(c(i-1, i, i, i-1), c(j-1, j-1, j, j))
			polygon(x[sub], y[sub], col=safesample(color),
				border=NA)
		}
	}
}


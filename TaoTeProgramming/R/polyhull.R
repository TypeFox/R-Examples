"polyhull" <-
function (num=100, lambda=20, size=.05, 
          color=grep("^gray", colors(), value=TRUE), seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	canvas()
	for(i in 1:num) {
		sides <- rpois(1, lambda) + 2
		x <- runif(1) + runif(sides, -size, size)
		y <- runif(1) + runif(sides, -size, size)
		ch <- chull(x, y)
		polygon(x[ch], y[ch], col=safesample(color), border=NA)
	}
}


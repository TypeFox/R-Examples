"butterflies" <-
function (num=c(100, 10), scale=c(.03, .1), 
	color=grep("^gray", colors(), value=TRUE), seed=NULL)
{
	if(length(seed)) set.seed(seed)

	n <- max(length(num), length(scale))
	num <- rep(num, length=n)
	scale <- rep(scale, length=n)

	canvas()
	for(j in 1:n) {
		sj <- scale[j]
		for(i in 1:num[j]) {
			butterfly(loc=runif(2, sj, 1-sj), scale=sj, 
				color=safesample(color))
		}
	}
}


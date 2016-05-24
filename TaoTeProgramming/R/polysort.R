"polysort" <-
function (num=100, lambda=20, size=.05, 
          color=grep("^gray", colors(), value=TRUE), seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	canvas()
	for(i in 1:num) {
		sides <- rpois(1, lambda) + 2
		polygon(runif(1) + sort(runif(sides, -size, size)), 
			runif(1) + sort(runif(sides, -size, size)), 
			col=safesample(color), border=NA)
	}
}


"polytile" <-
function (num=100, lambda=20, 
          color=grep("^gray", colors(), value=TRUE), seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	canvas()
	for(i in 1:num) {
		sides <- rpois(1, lambda) + 2
		polygon(runif(sides), runif(sides), 
			col=safesample(color), border=NA)
	}
}

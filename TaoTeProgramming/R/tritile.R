"tritile" <-
function (num=100, color=grep("^gray", colors(), value=TRUE), seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	canvas()
	for(i in 1:num) {
		polygon(runif(3), runif(3), col=safesample(color), border=NA)
	}
}


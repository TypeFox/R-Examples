"roads" <-
function (num=100, hprob=.6, color="black", lwd=1, seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	canvas()
	for(i in 1:num) {
		tri <- runif(3)
		if(runif(1) < hprob) {
			segments(tri[1], tri[2], tri[1], tri[3], 
				col=safesample(color), lwd=safesample(lwd))
		} else {
			segments(tri[1], tri[2], tri[3], tri[2], 
				col=safesample(color), lwd=safesample(lwd))
		}
	}
}


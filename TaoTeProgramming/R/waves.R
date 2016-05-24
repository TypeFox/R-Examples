"waves" <-
function (num=20, lwd=2, color=grep("^gray", colors(), value=TRUE), 
	seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	canvas()
	for(i in 1:num) for(j in 1:num) {
		ran <- runif(2)
		x <- seq(ran[1], ran[2], length=100)
		lines(x -.5 + i/num, .5 * tan(x * pi) + j/num,
			col=safesample(color), lwd=safesample(lwd))
	}
}


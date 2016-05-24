"bundle" <-
function (num=100, orig=NULL, width=.1, dest=NULL, 
          color="black", seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	if(!length(orig)) orig <- runif(2)
	if(!length(dest)) dest <- runif(2)

	width <- rep(width, length=4)

	segments(orig[1] + runif(num, -width[1], width[1]), orig[2] + 
		runif(num, -width[2], width[2]),
		dest[1] + runif(num, -width[3], width[3]), dest[2] + 
		runif(num, -width[4], width[4]), col=color)
}


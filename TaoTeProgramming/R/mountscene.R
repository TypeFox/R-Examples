"mountscene" <-
function (levels=5, vlen=200, tilt=.2, df=7, color=NULL, box=TRUE, seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	vlen <- rep(vlen, length=levels)
	df <- rep(df, length=levels)
	if(!length(color)) {
		color <- paste("gray", round(seq(1, 69, length=levels+2)[-c(1, 
			levels+2)]), sep="")
	} else {
		if(length(color) != levels) {
			stop("'color' should have length ", levels)
		}
	}
	lows <- seq(0, .7, length=levels)
	highs <- seq(.2, .9, length=levels)
	
	canvas()
	for(i in levels:1) {
		ttilt <- runif(1, -tilt, tilt)
		tiltvec <- seq(-ttilt, ttilt, length=vlen[i])
		mountain(runif(vlen[i], lows[i], highs[i]) + tiltvec, df=df[i],
			color=color[i])
	}
	if(box) polygon(c(0,0,1,1), c(0,1,1,0))
}


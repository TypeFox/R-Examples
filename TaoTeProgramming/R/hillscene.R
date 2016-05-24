"hillscene" <-
function (num=100, vlen=200, tilt=.2, df=7, color=NULL, seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	vlen <- rep(vlen, length=num)
	df <- rep(df, length=num)
	if(!length(color)) {
		color <- grep("^gray", colors(), value=TRUE)
	}
	lows <- seq(0, .7, length=num)
	highs <- seq(.2, .9, length=num)
	
	canvas()
	for(i in num:1) {
		ttilt <- runif(1, -tilt, tilt)
		tiltvec <- seq(-ttilt, ttilt, length=vlen[i])
		mountain(runif(vlen[i], lows[i], highs[i]) + tiltvec, df=df[i],
			color=safesample(color))
	}
	polygon(c(0,0,1,1), c(0,1,1,0))
}


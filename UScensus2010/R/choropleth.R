choropleth <-
function(sp,dem="P0010001",cuts=list("quantile",seq(0, 1, 0.25)),color=list(fun="hsv",attr=list(h = c(.4,.5,.6,.7), s = .6, v = .6, alpha=1)),main=NULL,sub="Quantiles (equal frequency)",border="transparent",legend=list(pos="bottomleft",title="Population Count"),type="plot",...){
	
m <- match.call()
m[[1]] <- as.name(paste("choropleth",type,sep="."))

eval(m, parent.frame())
}


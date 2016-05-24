dgumbel <-
function(x,scale=1,location=0,log=FALSE)
{
	fx <- 1/scale*exp(-(x-location)/scale-exp(-(x-location)/scale))
	if(log) return(log(fx))
	else return(fx)
}


dgp <-
function(x,shape=1,scale=1,log=FALSE)
{
	fx <- 1/scale*(1-shape*x/scale)^((1-shape)/shape)
	if(log) return(log(fx))
	else return(fx)
}


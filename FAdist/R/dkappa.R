dkappa <-
function(x,shape=1,scale=1,log=FALSE)
{
	fx <- shape/scale*(shape+(x/scale)^shape)^(-(shape+1)/shape)
	if(log) return(log(fx))
	else return(fx)
}


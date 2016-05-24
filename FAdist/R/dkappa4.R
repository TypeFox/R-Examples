dkappa4 <-
function(x,shape1,shape2,scale=1,location=0,log=FALSE)
{
	fx <- (1-shape1*(x-location)/scale)^(1/shape1-1)/scale*
		pkappa4(x,shape1,shape2,scale,location)^(1-shape2)
	if(log) 
		return(log(fx))
	else return(fx)
}
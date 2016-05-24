pkappa4 <-
function(q,shape1,shape2,scale=1,location=0,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- (1-shape2*(1-shape1/scale*(q-location))^(1/shape1))^(1/shape2)
	if(!lower.tail) 
	Fx <- 1 - Fx
	if(log.p) 
		Fx <- log(Fx)
	return(Fx)
}
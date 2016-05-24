qkappa4 <-
function(p,shape1,shape2,scale=1,location=0,lower.tail=TRUE,log.p=FALSE) 
{
	if(log.p) 
		p <- exp(p)
	if(!lower.tail) 
		p <- 1 - p
	xF <- location+scale/shape1*(1-((1-p^shape2)/shape2)^shape1)
	return(xF)
}

`pvfx` <-
function(test)
{
	if (test<0) test=0
	return(1-pchibar(test, df=c(0,1), wt=c(.5,.5)))
}


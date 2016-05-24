`scaleBranchingtimes` <-
function(x, basal = 100)
{
	x <- rev(sort(as.numeric(x)));
	sf <- basal/max(x);
	return(sf*x);

}


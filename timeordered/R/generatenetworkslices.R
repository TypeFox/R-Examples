generatenetworkslices <-
function(g, timedeltas)
{

	
	return(apply(timedeltas, 1, function(x) { generatetimeaggregatednetwork(g, x[1],x[2]) }))
}


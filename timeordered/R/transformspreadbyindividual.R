transformspreadbyindividual <-
function(sa)
{
	# returns mean by individual at each timepoint
	return(sapply(by(sa, sa$startvertex, function(x) { apply(x[, -1, drop=FALSE], 2, mean) }), function(x){x}))
}


shortesthoppath <-
function(g, startvertexname, startvertextime, stopvertexname, stopvertextime)
{
	if (length(startvertexname) != 1 | length(stopvertexname) != 1)
	{
		stop("must provide single startvertex and stop vertex")	
	}
	startvertex <- V(g)[V(g)$Name==startvertexname & V(g)$Time==startvertextime]
	stopvertex <- V(g)[V(g)$Name==stopvertexname & V(g)$Time==stopvertextime]
	vertices <- get.shortest.paths(g, startvertex, stopvertex, mode="out", weights=E(g)$HopCost)
        if (is.list(vertices) && "vpath" %in% names(vertices)) { vertices <- vertices$vpath }
	
	return(V(g)[vertices[[1]] ])
}


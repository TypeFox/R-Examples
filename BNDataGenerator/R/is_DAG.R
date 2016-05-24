is_DAG = function (arcs_mat) 
{
	unmakeMG = function (arcs_mat) 
	{
		d = nrow(arcs_mat)
		ug = dg = bg = arcs_mat
		M = expand.grid(dg = 0:1, ug = 0:1, bg = 0:1)
		i = strtoi(as.character(arcs_mat), 2)
		GG = M[i + 1, ]
		ug[, ] = GG[, 2]
		dg[, ] = GG[, 1]
		bg[, ] = GG[, 3]
		if (any(ug != t(ug))) 
			stop("Undirected edges are wrongly coded.")
		if (any(bg != t(bg))) 
			stop("Undirected edges are wrongly coded.")
		return(list(dg = dg, ug = ug, bg = bg))
	}

    comp = unmakeMG(arcs_mat)
    ug = comp$ug
    dag = comp$dg
    bg = comp$bg
    
	out = TRUE
	
    if (any(arcs_mat > 100)) {
        warning("There are double edges.")
        out = FALSE
    }
	
    if (!is_acyclic(dag)) {
        warning("Not acyclic.")
        out = FALSE
    }
    return(out)
}
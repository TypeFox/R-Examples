fromto_to_mat = function(fromto, nodename)
{
	if(dim(fromto)[1] == 0)
	{
		stop("It has not any arc");
	}
	
	num_of_nodes = length(nodename)
	arcs_mat = matrix(0, num_of_nodes, num_of_nodes)

	arcs_order_mat = cbind(nodename, c(1:length(nodename)))
	temp_arcs = cbind(match(fromto[,1], arcs_order_mat), match(fromto[,2], arcs_order_mat))


	if (length(temp_arcs) > 0)
	{
		for (i in 1:dim(temp_arcs)[1])
		{
			from = as.numeric(temp_arcs[i,1])
			to = as.numeric(temp_arcs[i,2])
			arcs_mat[from, to] = arcs_mat[from, to] + 1
		}
	}
	
	dimnames(arcs_mat)[[1]] = nodename
	dimnames(arcs_mat)[[2]] = nodename
	
	return(arcs_mat);
}
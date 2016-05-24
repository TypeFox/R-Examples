# Written by Jae-seong Yoo 20141222

make_topology = function (nodes, topology, Probs=NULL, nodename=NULL, cardinality=NULL)
{
	# Check Num of Nodes
	NeedMoreNodes = function(num_of_nodes)
	{
		if(nodes < num_of_nodes)
			stop("Need More Nodes!");
	}
	
	switch (topology,
		"Collapse" = {NeedMoreNodes(3)},
		"Line" = {NeedMoreNodes(3)},
		"Star" = {NeedMoreNodes(3)},
		"PseudoLoop" = {NeedMoreNodes(3)},
		"Diamond" = {NeedMoreNodes(4)},
		"Rhombus" = {NeedMoreNodes(4)},
	)
	
	
	# arcs_mat_mat
	arcs_mat = matrix(0, nodes, nodes);
	switch (topology,
		"Collapse" = {	arcs_mat[,nodes] = 1;
								arcs_mat[nodes,nodes] = 0;
							},
		"Line" = {	for (i in 1:(nodes-1))
						{
							arcs_mat[i,(i+1)] = 1
						}
					},
		"Star" = {	arcs_mat[1,] = 1
						arcs_mat[1,1] = 0
					},
		"PseudoLoop" = {	arcs_mat[1, nodes] = 1
									for (i in 1:(nodes-1))
									{
										arcs_mat[i,(i+1)] = 1
									}
								},
		"Diamond" = {	arcs_mat[1,] = 1
								arcs_mat[1,1] = 0
								arcs_mat[,nodes] = 1
								arcs_mat[1,nodes] = 0
								arcs_mat[nodes,nodes] = 0
							},
		"Rhombus" = {	arcs_mat[1,] = 1
								arcs_mat[2,] = 1
								arcs_mat[(1:2),(1:2)] = 0
								arcs_mat[nodes,nodes] = 0
							},
	)
	
	
	# Check Input Probs & cardinality
	checker = check_cardinality(arcs_mat = arcs_mat, nodename = nodename, cardinality = cardinality)
	cardinality = checker$cardinality;
	num_of_probs = checker$num_of_probs;
	nodename = checker$nodename;
	
	
	# Random Probs
	if (is.null(Probs) & is.null(cardinality))
	{
		Probs = list()
		switch (topology,
			"Collapse" = {	for (i in 1:(nodes-1))
									{
										Probs[[i]] = runif(1)
									}
									Probs[[nodes]] = runif(2^(nodes-1))
								},
			"Line" = {	Probs[[1]] = runif(1)
							for (i in 2:nodes)
							{
								Probs[[i]] = runif(2)
							}
						},
			"Star" = {	Probs[[1]] = runif(1)
							for (i in 2:nodes)
							{
								Probs[[i]] = runif(2)
							}
						},
			"PseudoLoop" = {	Probs[[1]] = runif(1)
										for (i in 2:(nodes-1))
										{
											Probs[[i]] = runif(2)
										}
										Probs[[nodes]] = runif(4)
									},
			"Diamond" = {	Probs[[1]] = runif(1)
									for (i in 2:(nodes-1))
									{
										Probs[[i]] = runif(2)
									}
									Probs[[nodes]] = runif(2^(nodes-2))
								},
			"Rhombus" = {	Probs[[1]] = runif(1)
									Probs[[2]] = runif(1)
									for (i in 3:nodes)
									{
										Probs[[i]] = runif(2^2)
									}
								},
		)
	} else if (is.null(Probs)) {
		Probs = list()
		for (i in 1:length(num_of_probs))
		{
			Probs[[i]] = runif(num_of_probs[i])
		}
	}
	
	
	result = list(	arcs_mat = arcs_mat,
						Probs = Probs,
						nodename = nodename,
						cardinality = cardinality,
						num_of_nodes = nodes
					)
	return(result)
}
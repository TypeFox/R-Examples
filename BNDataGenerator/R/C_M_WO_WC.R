C_M_WO_WC = function(target_arcs_mat, learnt_arcs_mat)
{
	nodes = dim(target_arcs_mat)[2]
	
	C = 0;	M = 0;	WO = 0;	WC = 0;

	for (i in 1:nodes)
	{
		C = C + abs(			sum((target_arcs_mat[,i] == 1) & (learnt_arcs_mat[,i] == 1)
											& (target_arcs_mat[,i] == learnt_arcs_mat[,i])))
		M = M + abs(		sum((target_arcs_mat[,i] == 1) & (learnt_arcs_mat[,i] == 0)
											& (target_arcs_mat[,i] != learnt_arcs_mat[,i])) -
									sum((target_arcs_mat[,i] == 1) & (learnt_arcs_mat[i,] == 1)
											& (target_arcs_mat[,i] == learnt_arcs_mat[i,])))
		WO = WO + abs(	sum((target_arcs_mat[,i] == 1) & (learnt_arcs_mat[i,] == 1)
											& (target_arcs_mat[,i] == learnt_arcs_mat[i,])))
		WC = WC + abs(	sum((target_arcs_mat[,i] == 0) & (learnt_arcs_mat[,i] == 1)
											& (target_arcs_mat[,i] != learnt_arcs_mat[,i])) -
									sum((target_arcs_mat[i,] == 0) & (learnt_arcs_mat[i,] == 1)
											& (target_arcs_mat[i,] != learnt_arcs_mat[i,])))
	}

	result = t(as.matrix(c(C, M, WO, WC)))
	dimnames(result)[[2]] = c("C", "M", "WO", "WC")

	return(result)
}

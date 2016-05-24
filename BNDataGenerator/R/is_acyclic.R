is_acyclic = function (arcs_mat) 
{
	
	transClos = function (arcs_mat) 
	{
		if (nrow(arcs_mat) == 1) 
			return(arcs_mat)
		A = arcs_mat
		diag(A) = 1
		repeat {
			B = sign(A %*% A)
			if (all(B == A)) 
				break
			else A = B
		}
		diag(A) = 0
		A
	}

	B = transClos(arcs_mat)
	l = B[lower.tri(B)]
	u = t(B)[lower.tri(t(B))]
	com = (l & u)
	return(all(!com))
}
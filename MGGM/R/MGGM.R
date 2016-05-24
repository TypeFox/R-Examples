MGGM.path <- function(S_bar, # p by p*L matrix: [S_1, ... ,S_L] with S_l being the l-th sample cov matrix
			nn, # L dimensional vector storing the sample sizes
			Lambda1.vec, Lambda2.vec, # lambda grids for lambda1 and lambda2
			graph, # matrix specifying pairs of precision matrices to be penalized to be similar; see example.R for constructions 
			tau = .01, MAX_iter=200, eps_mat = 1e-4){
	p = dim(S_bar)[1]
	L = dim(S_bar)[2] / p
	grid.lambda1 = length(Lambda1.vec)
	grid.lambda2 = length(Lambda2.vec)
	# four matrices for input #
	covmat_inverse_path = matrix(rep(diag(p),L),p,p*L*grid.lambda1*grid.lambda2)
	covmat_inverse_con_path = covmat_inverse_path
	covmat_path = covmat_inverse_path
	covmat_con_path = covmat_inverse_path
	NumOfEdge = dim(graph)[2]
	
	out <- .C("matrix_grouping_path", S_bar = as.double(S_bar),
		covmat_inverse_path=as.double(covmat_inverse_path),
	  	covmat_path = as.double(covmat_path),
		covmat_inverse_con_path=as.double(covmat_inverse_con_path),
	  	covmat_con_path = as.double(covmat_con_path),
	  	Lambda1=as.double(Lambda1.vec),Lambda2 = as.double(Lambda2.vec),
	  	Tau=as.double(tau), grid_lambda1 = as.integer(grid.lambda1), 
	  	grid_lambda2 = as.integer(grid.lambda2), 
	  	Graph=as.integer(graph),sample_size=as.double(nn),
	  	pp=as.integer(p),LL=as.integer(L),
	  	NumOfEdges=as.integer(NumOfEdge),MAX_DC_ITER0=as.double(MAX_iter),
		  eps_mat=as.double(eps_mat),PACKAGE = "MGGM")
	
	sol_path = list()
	sol_path$sol_nonconvex = array(out$covmat_inverse_path, dim=c(p,p*L,grid.lambda2,grid.lambda1))
	sol_path$sol_convex = array(out$covmat_inverse_con_path, dim=c(p,p*L,grid.lambda2,grid.lambda1))
	return (sol_path)
}

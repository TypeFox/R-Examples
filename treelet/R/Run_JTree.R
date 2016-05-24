Run_JTree <-
function(X, maxlev,whichsave){

if(maxlev>=ncol(X)){

	message("can only compute up to p-1 merges")

	} else{

    message("computing the Correlation....")
	C = X
	cc = cov2cor(C)
    message("building the tree......")
    maxTree = Build_JTree(C, cc, maxlev,whichsave)
    message("computing the basis for the level......")
    out = JTree_Basis(maxTree$Zpos, maxTree$T, maxTree$PCidx, 
        maxlev, maxTree$all_nodes,whichsave)

        return(list(basis = out$basis, Zpos = maxTree$Zpos, 
        T = maxTree$T, PCidx = maxTree$PCidx, 
        all_nodes = maxTree$all_nodes, TreeCovs=maxTree$TreeCovs))

	}
}

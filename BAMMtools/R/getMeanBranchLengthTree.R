#############################################################
#
#	getMeanBranchLengthTree(....)
#
 
getMeanBranchLengthTree <- function(ephy, rate='speciation') {
	
	if (class(ephy) == 'bammdata') {
		v <- as.phylo.bammdata(ephy);
	}
	obj <- getMarginalBranchRateMatrix(ephy,verbose=FALSE);

	if (ephy$type == 'diversification'){
		if (rate == 'speciation'){
			el <- rowMeans(obj$lambda_branch_matrix);
		}else if (rate == 'extinction'){
			el <- rowMeans(obj$mu_branch_matrix);
		}else if (rate == 'ndr'){
			el  <- rowMeans(obj$lambda_branch_matrix) - rowMeans(obj$mu_branch_matrix);				
		}else{
			stop("invalid rate specification in getMeanBranchLengthTree");
		}
		
	}else if (ephy$type == 'trait'){
		el <- rowMeans(obj$beta_branch_matrix);
	
	}else{
		stop("error in getMeanBranchLengthTree - \nproblem with supplied ephy object");
	}
	
	v$edge.length <- el;
	tmp <- list();
	tmp$phy <- v;
	tmp$median <- median(el);
	tmp$mean <- mean(el);		
		
	
	return(tmp);		
}

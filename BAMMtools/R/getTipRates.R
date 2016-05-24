#############################################################
#
#	getTipRates(....)
#
# Returns a list with:
# lambda = matrix of tip rates where rows are species and columns are posterior samples, 
# mu if ephy$type == 'diversification',
# beta if ephy$type='trait',
# lambda.avg, mu.avg, beta.avg: named vector of average tip rates. 
# If returnNetDiv = TRUE, then a matrix and average vector for net div rates is returned.

getTipRates <- function(ephy, returnNetDiv = FALSE, statistic = 'mean') {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}
	if (!statistic %in% c('mean','median')) {
		stop("statistic must be either 'mean' or 'median'.");
	}

	obj <- list();
	if (ephy$type == 'diversification') {
		if (returnNetDiv) {
			obj$netdiv <- do.call(cbind, ephy$tipLambda) - do.call(cbind, ephy$tipMu);
			rownames(obj$netdiv) <- as.phylo.bammdata(ephy)$tip.label;
			
			if (statistic == 'mean') {
				obj$netdiv.avg <- rowMeans(obj$netdiv);
			}
			if (statistic == 'median') {
				obj$netdiv.avg <- apply(obj$netdiv, 1, median);
			}
		}
		
		if (!returnNetDiv) {
			obj$lambda <- do.call(cbind, ephy$tipLambda);
			rownames(obj$lambda) <- as.phylo.bammdata(ephy)$tip.label;
		
			obj$mu <- do.call(cbind, ephy$tipMu);
			rownames(obj$mu) <- as.phylo.bammdata(ephy)$tip.label;
		
			if (statistic == 'mean') {
				obj$lambda.avg <- rowMeans(obj$lambda);
				obj$mu.avg <- rowMeans(obj$mu);
			}
			if (statistic == 'median') {
				obj$lambda.avg <- apply(obj$lambda, 1, median);
				obj$mu.avg <- apply(obj$mu, 1, median);
			}
		}
	}
	
	if (ephy$type == 'trait') {
		obj$beta <- do.call(cbind, ephy$tipLambda);
		rownames(obj$beta) <- as.phylo.bammdata(ephy)$tip.label;
		
		if (statistic == 'mean') {
			obj$beta.avg <- rowMeans(obj$beta);
		}
		if (statistic == 'median') {
			obj$betabeta.avg <- apply(obj$beta, 1, median);
		}
	}
	return(obj);
}










#  File networksis/R/simulate.sisnetwork.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
#########################################################################

simulate.sisnetwork <- function(object, nsim = 1, seed = NULL, save.networks = FALSE, ...)
{
	if(class(object) == "network")
	{
	    object <- network::as.sociomatrix(object)
    	object <- list(rows = apply(object, 1, sum), cols = apply(object, 2, sum))
	    class(object) <- "sisnetwork"
	} else if(!("sisnetwork" %in% class(object)))
	{
		stop("simulate.sisnetwork requires an object of class 'sisnetwork' or 'network'.")
	}

	if(save.networks)
		out.list <- vector("list", nsim)

	# Set seed if specified
	if(!is.null(seed))
		set.seed(seed)

	# Initialize the vector of probabilities
	probvec <- rep(0, nsim)

	# Extract row sums and column sums, and rearrange column sums in descending order
	row.sums <- object$rows
	col.sums <- object$cols
	order.col.sums <- order(-col.sums)
	col.sums <- col.sums[order.col.sums]

	# Determine dimensions of matrix and initialize matrix for new graph(s)
	nrow <- length(row.sums)
	ncol <- length(col.sums)
	mat <- matrix(0, nrow = nrow, ncol = ncol)
	
	# Calculate number of edges in graph
	nedges <- sum(row.sums)

	nsim.in <- 1
	
	# Keep only non-zero row and column sums.
	col.sums <- col.sums[col.sums > 0]
	row.sums <- row.sums[row.sums > 0]
	nrow.nonnull <- length(row.sums)
	ncol.nonnull <- length(col.sums)

	if(ncol.nonnull == 0 & nrow.nonnull == 0)
		stop("This is an empty graph.  simulate.sisnetwork requires at least one edge.")
		
	for(i in 1 : nsim)
	{
		s <- .C("sissim", 
				rowsums = as.integer(row.sums),
				colsums = as.integer(col.sums),
				newmat = integer(nrow.nonnull * ncol.nonnull), 
				nrow = as.integer(nrow.nonnull), 
				ncol = as.integer(ncol.nonnull), 
				as.double(nsim.in),
				prob = as.double(1), 
				probvec = double(nsim.in),
				PACKAGE = "networksis")
   
		probvec[i] <- s$probvec
		if(save.networks)
		{
			mat[1 : nrow.nonnull, 1 : ncol.nonnull] <- matrix(s$newmat, ncol = ncol.nonnull)
			mat <- mat[, order(order.col.sums)]
			if(nsim > 1)
			{
				out.list[[i]] <- network::as.network.matrix(mat, matrix.type = "bipartite")
			} else{
				out.list <- network::as.network.matrix(mat, matrix.type = "bipartite")
			}
		}
	}

	log.graphspace.mean <- log(mean(1 / exp(probvec)))
	log.graphspace.SE.mean <- log(var(1 / exp(probvec)) ^ .5 / nsim ^ .5)
	meanprob <- mean(-probvec)
	varprob <- var(probvec)
	g <- function(s)
	{
		k <- 1 : 100
		fk <- k - k + nsim / (nsim + 1)
		
		for(kk in k[-length(k)])
		{
			fk[(kk+1) : max(k)] <- nsim * nsim * fk[kk : (max(k) - 1)] / ((kk + 1) * (nsim + 1) * (nsim + 2 * kk))
		}
		1 + sum((s ^ k) * fk)
	}
	log.graphspace.size.lne <- meanprob + log(g(0.5 * varprob))
	log.graphspace.SE.lne <- 0.5 * (2 * meanprob + varprob + log(varprob + 0.5 * varprob * varprob * (1 + (varprob + varprob * varprob / 4) / nsim)) - log(nsim))

	if(save.networks)
	{
		out.list <- list(networks = out.list, 
						log.prob = probvec,
						log.graphspace.size = log.graphspace.mean, 
						log.graphspace.SE = log.graphspace.SE.mean,
						log.graphspace.size.lne = log.graphspace.size.lne,
						log.graphspace.SE.lne = log.graphspace.SE.lne) 
		class(out.list) <- "network.series"
	}else{
		mat[1 : nrow.nonnull, 1 : ncol.nonnull] <- matrix(s$newmat, ncol = ncol.nonnull)[, order(order.col.sums)]
		out.list <- network::as.network.matrix(mat, matrix.type = "bipartite")
		out.list <- list(networks = out.list, 
						log.prob = probvec,
						log.graphspace.size = log.graphspace.mean, 
						log.graphspace.SE = log.graphspace.SE.mean,
						log.graphspace.size.lne = log.graphspace.size.lne,
						log.graphspace.SE.lne = log.graphspace.SE.lne) 
		class(out.list) <- "network.series"
	}

	return(out.list)
}

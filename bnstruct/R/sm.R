# test.compute.counts.nas <- function( data, bins, nas )
# {
# 	.Call("compute_counts_nas", data, bins, nas, package = "bnstruct")
# }
# 
# test.na.rows.int <- function( data )
# {
# 	.Call("na_rows_int", data, package = "bnstruct")
# }

# Silander and Myllymaki complete search, BDeu score, high memory occupation.
sm <- function(x, node.sizes, scoring.func = 0, cont.nodes = NULL, max.fanin = NULL, 
	layering = NULL, max.fanin.layers = NULL, ess = 1, cpc.mat = NULL ) 
{
	n.cases <- nrow(x)
	n.nodes <- ncol(x)
	
	storage.mode(node.sizes) <- "integer" # just to be sure
  storage.mode(scoring.func) <- "integer"
	
	# if no max.fanin is given, assume maximum possible
	if( is.null(max.fanin) )
		max.fanin <- n.nodes - 1
	
	# if no layering is given, assume one layer
	if( is.null(layering) )
	{
		layering <- rep( 1, n.nodes )
		n.layers <- 1
		if( !is.null(max.fanin.layers) )
			stop("no max.fanin.layers without layering.")
		else
			max.fanin.layers <- array(max.fanin)
	}
	else
	{
		n.layers <- length(unique(layering))
		if( is.null(max.fanin.layers) )
		{
      # derive max.fanin.layers from layering
		  max.fanin.layers <- array(cumsum(table(layering)))
      max.fanin.layers[1] <- 0 # default, no parents for nodes at layer 1
		}
	}

  if( ( !is.array(max.fanin.layers) ) || 
        ( dim(max.fanin.layers)[1] != n.layers ) || 
        ( is.matrix(max.fanin.layers) && dim(max.fanin.layers)[2] != n.layers) )
	  stop("max.fanin.layers must be either an array of #layers elements or a #layers x #layers squared matrix")
	
	# get matrix from 1D array
	if( !is.matrix(max.fanin.layers) && n.layers > 1 )
	{
		max.fanin.layers <- diag(max.fanin.layers)
		max.fanin.layers[upper.tri(max.fanin.layers)] <- -1 # don't care
	}	
	
	# thresholding with max.fanin
	max.fanin.layers[max.fanin.layers > max.fanin] <- max.fanin
	
	# trick, to avoid stupid R behaviour
	dim(max.fanin.layers) <- c(n.layers,n.layers)
	  
	# quantize with NAs
	levels <- rep( 0, n.nodes )
	levels[cont.nodes] <- node.sizes[cont.nodes]
	# data <- quantize.with.na.matrix( x, levels )
	data <- quantize.matrix( x, levels )	
  
	ifm <- impossible.family.mask( n.nodes, layering, max.fanin.layers )
  
  # remove parents not in cpc, if cpc.matrix is given
	if( !is.null(cpc.mat) )
    for( i in 1:n.nodes )
      ifm[ i, (.Call("fumt_mask", n_elements = n.nodes, pattern = which(cpc.mat[i,]==0), 
                     PACKAGE = "bnstruct") > 0) ] <- FALSE
  
  # aflml <- all.families.log.marginal.likelihood( data, node.sizes, ifm, ess )
	aflml <- .Call("all_fam_log_marg_lik", data, node.sizes, ifm, ess, scoring.func, PACKAGE = "bnstruct" )
	
	# bps <- find.best.parents( aflml )
	bps <- .Call("fbp", aflml = aflml, PACKAGE = "bnstruct");
	
	# sinks <- find.best.sinks( bps, aflml )
	sinks <- .Call("fbs", bps = bps, aflml = aflml, PACKAGE = "bnstruct");
	
	order <- find.best.ordering( sinks )
	
	parents <- find.best.network( order, bps )
	
	dag <- parents.to.dag( parents )
	
	return( dag )
}#SM

# quantize.with.na <- function(data, levels) 
# {
# 	nr <- nrow(data)
# 	nc <- ncol(data)
# 	
# 	quant <- data
# 	
# 	for( i in 1:nc )
# 	{
# 		if( levels[i] > 0 )	# to be discretized
# 		{
# 			quantiles <- quantile( data[,i], probs = (0:levels[i])/levels[i], na.rm = TRUE )
# 			# cut the range using the quantiles as break points.
# 			quant[,i] <- as.data.frame( cut( data[,i], quantiles, labels=FALSE, include.lowest=TRUE) )
# 		}
# 	}
# 	
# 	return(quant)
# }

# quantize.with.zeros.matrix <- function(data, levels) 
# {
# 	nr <- nrow(data)
# 	nc <- ncol(data)
# 	
# 	quant <- matrix(0,nr,nc)
# 	
# 	for( i in 1:nc )
# 	{
# 		if( levels[i] == 0 )	#already discrete
# 			quant[,i] <- as.matrix(data[,i],nr,1)
# 		else
# 		{
# 			quantiles <- quantile( data[,i], probs = (0:levels[i])/levels[i], na.rm = TRUE )
# 			# cut the range using the quantiles as break points.
# 			quant[,i] <- as.matrix( cut( data[,i], quantiles, labels=FALSE, include.lowest=TRUE),nr,1 )
# 		}
# 	}
# 	
# 	# turn NAs to 0s
# 	quant[is.na(quant)] <- 0
# 	return(quant)
# }

impossible.family.mask <- function( n.nodes, layering, max.fanin.layers)
{
	# filter based on maxFanIn
	# base.mask <- rep(0,2^n.nodes)
	# base.mask[ 2^((1:n.nodes)-1)+1 ] <- 1
	# base.mask <- fumt(base.mask) > max(diag(max.fanin.layers))
	base.mask <- (.Call("fumt_mask", n_elements = n.nodes, pattern = seq_len(n.nodes),
		PACKAGE = "bnstruct") > max(diag(max.fanin.layers)) );
	
	ifm <- !matrix( rep(base.mask,n.nodes), n.nodes, 2^n.nodes, byrow = TRUE)
	
	for( i in 1:n.nodes )
	{
		# mask <- rep(0, 2^n.nodes)
		# mask[ 2^(i-1)+1 ] <- 1
		# mask <- fumt(mask) > 0
		mask <- (.Call("fumt_mask", n_elements = n.nodes, 
			pattern = i, PACKAGE = "bnstruct") > 0);
		
		ifm[ i, mask ] <- FALSE
	}
	
	if( length(layering) > 1 )
	{
		for( i in 1:n.nodes )
		{
			# false entries if they violate layering
			invalidParents <- which( layering > layering[i] )
			# mask <- rep(0, 2^n.nodes)
			# mask[ 2^(invalidParents-1)+1 ] <- 1
			# mask <- fumt(mask) > 0
			mask <- (.Call("fumt_mask", n_elements = n.nodes, 
				pattern = invalidParents, PACKAGE = "bnstruct") > 0);
			
			ifm[ i, mask ] <- FALSE
			
			# false entries if they violate max.fanin.layers-layering
			# this refines the filter on max.fanin above 
			# (which was not specific to a particular child)
			for ( j in unique(layering[layering <= layering[i]]) )
			{
				if( max.fanin.layers[j,layering[i]]==-1 ) 
					next # no restriction

				if( j == layering[i] )
					parents <- (layering<=j)
				else
					parents <- (layering==j)

				parents[i] <- FALSE
				parents <- which(parents)
				# mask <- rep(0, 2^n.nodes)
				# mask[ 2^(parents-1)+1 ] <- 1
				# mask <- fumt(mask) > max.fanin.layers[j, layering[i]];
				mask <- ( .Call("fumt_mask", n_elements = n.nodes, pattern = parents,
					PACKAGE = "bnstruct") > max.fanin.layers[j, layering[i]] );
				ifm[i, mask] <- FALSE;
			}
		}
	}
	
	return(ifm)
}

# Fast Upwards Moebius Transform 
# Based on code by Smets, http://iridia.ulb.ac.be/~psmets/
# fumt <- function( h0 )
# {
# 	lh0 <- length( h0 )
# 	n <- round( log2( lh0 ) )
# 	# if( lh0 != 2^n)
# 	#	stop("length of h0 must be a power of 2")
# 	
# 	hn <- h0
# 	
# 	for( i in seq_len(n) )
# 	{
# 		dim(hn) <- c( 2^(i-1), 2^(n+1-i) )
# 		ind <- seq_len( 2^(n-i) ) * 2
# 		hn[,ind] <- hn[,ind] + hn[,ind-1] 
# 	}
# 	
# 	dim(hn) <- lh0
# 	return( hn )
# }

# all.families.log.marginal.likelihood <- function( data, node.sizes, ifm, ess )
# {
# 	n.nodes <- ncol(data)
# 	LM <- matrix( -Inf, n.nodes, 2^n.nodes )
# 	bitmask <- 2^(0:(n.nodes-1))
# 	
# 	for( i in 1:n.nodes )
# 	{
# 		# select possible parent sets and compute log-likelihood
# 		possible.families <- which(ifm[i,])
# 		for( k in possible.families )
# 		{
# 			pa <- which( bitAnd(k-1,bitmask) > 0 )
# 			# LM[i,k] <- log.likelihood.na( pa, i, node.sizes, ess, data)
# 			LM[i,k] <- log.lik.na( node.sizes[c(i,pa)], ess, data[,c(i,pa)] ) 
# 		}
# 	}
# 	
# 	return(LM)
# }

# log.lik.na <- function( node.sizes, ess, data )
# {
# 	n.nodes <- length(node.sizes)
# 	
# 	if( n.nodes > 1 )
# 	{
# 		prod.sizes.pa <- prod(node.sizes[2:n.nodes])
# 		prod.sizes <- prod.sizes.pa * node.sizes[1]
# 		na.rows <- .Call("na_rows_int", mat = data, PACKAGE = "bnstruct")
# 		
# 		# counts <- compute.counts( data, node.sizes )
# 		# prior <- array( ess/prod.sizes, node.sizes ); # faster !
# 		counts <- .Call("compute_counts_nas", data, node.sizes, na.rows, package = "bnstruct")
# 		prior <- ess/prod.sizes;
# 		
# 		# correct for NAs with maximum a posteriori estimate
# 		n.na <- sum(na.rows)
# 		counts <- counts + n.na * (counts + prior) / ( nrow(data) - n.na + ess )
# 		
# 		# LL <- sum( lgamma(prior+counts) - lgamma(prior) )
# 		LL <- sum(lgamma(prior+counts)) - prod.sizes*lgamma(prior) # faster!
# 		# alpha_ij <- colSums( prior )
# 		# alpha_ij <- array( ess/prod.sizes.pa, node.sizes[2:n.nodes] )
# 		alpha_ij <- ess/prod.sizes.pa
# 		N_ij <- colSums( counts )
# 		# return( LL + sum( lgamma(alpha_ij) - lgamma(alpha_ij+N_ij) ) )
# 		return( LL + prod.sizes.pa*lgamma(alpha_ij) - sum(lgamma(alpha_ij+N_ij)) ) # faster !
# 	}
# 	else
# 	{
# 		na.rows <- as.integer(is.na(data))
# 		# counts <- compute.counts( data, node.sizes )
# 		counts <- .Call("compute_counts_nas", data, node.sizes, na.rows, package = "bnstruct")
# 		#prior <- rep( ess/node.sizes, node.sizes );
# 		prior <- ess/node.sizes
# 
# 		# correct for NAs with maximum a posteriori estimate
# 		n.na <- sum(na.rows) 
# 		counts <- counts + n.na * (counts + prior) / ( length(data) - n.na + ess )
# 
# 		# LL <- sum( lgamma(prior+counts) - lgamma(prior) )
# 		LL <- sum(lgamma(prior+counts)) - node.sizes*lgamma(prior) # faster!
# 		return( LL - lgamma(1+sum(counts)) )
# 	}
# }

log.lik <- function( node.sizes, ess, data )
{
  n.nodes <- length(node.sizes)
  
  if( n.nodes > 1 )
  {
    prod.sizes.pa <- prod(node.sizes[2:n.nodes])
    prod.sizes <- prod.sizes.pa * node.sizes[1]
    # na.rows <- .Call("na_rows_int", mat = data, PACKAGE = "bnstruct")
    
    # counts <- compute.counts( data, node.sizes )
    # prior <- array( ess/prod.sizes, node.sizes ); # faster !
    # counts <- .Call("compute_counts_nas", data, node.sizes, na.rows, package = "bnstruct")
    counts <- .Call("compute_counts", data, node.sizes, PACKAGE = "bnstruct")
    prior <- ess/prod.sizes;
    
    # correct for NAs with maximum a posteriori estimate
    # n.na <- sum(na.rows)
    # counts <- counts + n.na * (counts + prior) / ( nrow(data) - n.na + ess )
    
    # LL <- sum( lgamma(prior+counts) - lgamma(prior) )
    LL <- sum(lgamma(prior+counts)) - prod.sizes*lgamma(prior) # faster!
    # alpha_ij <- colSums( prior )
    # alpha_ij <- array( ess/prod.sizes.pa, node.sizes[2:n.nodes] )
    alpha_ij <- ess/prod.sizes.pa
    N_ij <- colSums( counts )
    # return( LL + sum( lgamma(alpha_ij) - lgamma(alpha_ij+N_ij) ) )
    return( LL + prod.sizes.pa*lgamma(alpha_ij) - sum(lgamma(alpha_ij+N_ij)) ) # faster !
  }
  else
  {
    # na.rows <- as.integer(is.na(data))
    # counts <- compute.counts( data, node.sizes )
    # counts <- .Call("compute_counts_nas", data, node.sizes, na.rows, package = "bnstruct")
    # prior <- rep( ess/node.sizes, node.sizes );
    counts <- .Call("compute_counts", data, node.sizes, PACKAGE = "bnstruct")
    prior <- ess/node.sizes
    
    # correct for NAs with maximum a posteriori estimate
    # n.na <- sum(na.rows) 
    # counts <- counts + n.na * (counts + prior) / ( length(data) - n.na + ess )
    
    # LL <- sum( lgamma(prior+counts) - lgamma(prior) )
    LL <- sum(lgamma(prior+counts)) - node.sizes*lgamma(prior) # faster!
    return( LL - lgamma(1+sum(counts)) )
  }
}


# log.likelihood.na <- function( pa, ni, node.sizes, ess, data )
# {
# 	n.pa <- length(pa)
# 	prod.sizes.pa <- prod(node.sizes[pa])
# 	prod.sizes <- prod.sizes.pa * node.sizes[ni]
# 	
# 	counts <- compute.counts( data[,c(ni,pa)], node.sizes[c(ni,pa)] )
# 	prior <- array( ess/prod.sizes, node.sizes[c(ni,pa)] );
# 	if( n.pa > 0 )
# 	{
# 		# correct for NAs with maximum a posteriori estimate
# 		
# 		n.na <- sum(.Call("na_rows", mat =data[,c(pa,ni)], PACKAGE = "bnstruct"))
# 		counts <- counts + n.na * (counts + prior) / ( nrow(data) - n.na + ess )
# 	
# 		# LL <- sum( lgamma(prior+counts) - lgamma(prior) )
# 		LL <- sum(lgamma(prior+counts)) - prod.sizes*lgamma(ess/prod.sizes) # faster!
# 		# alpha_ij <- colSums( prior )
# 		alpha_ij <- array( ess/prod.sizes.pa, node.sizes[pa] );
# 		N_ij <- colSums( counts )
# 		# return( LL + sum( lgamma(alpha_ij) - lgamma(alpha_ij+N_ij) ) )
# 		return( LL + prod.sizes.pa*lgamma(ess/prod.sizes.pa) - sum(lgamma(alpha_ij+N_ij)) ) # faster !
# 	}
# 	else
# 	{
# 		# correct for NAs with maximum a posteriori estimate
# 		n.na <- sum(is.na(data[,ni])) 
# 		counts <- counts + n.na * (counts + prior) / ( nrow(data) - n.na + ess )
# 	
# 		# LL <- sum( lgamma(prior+counts) - lgamma(prior) )
# 		LL <- sum(lgamma(prior+counts)) - prod.sizes*lgamma(ess/prod.sizes) # faster!
# 		return( LL - lgamma(1+sum(counts)) )
# 	}
# }

# log.likelihood.na.rows <- function( pa, ni, node.sizes, ess, data )
# {
# 	n.pa <- length(pa)
# 	prod.sizes.pa <- prod(node.sizes[pa])
# 	prod.sizes <- prod.sizes.pa * node.sizes[ni]
# 	
# 	counts <- compute.counts( data[,c(pa,ni)], node.sizes[c(pa,ni)] )
# 	prior <- array( ess/prod.sizes, node.sizes[c(pa,ni)] );
# 	if( n.pa > 0 )
# 	{
# 		# correct for NAs with maximum a posteriori estimate
# 		n.na <- sum(rowSums(is.na(data[,c(pa,ni)]))>0)
# 		counts <- counts + n.na * (counts + prior) / ( nrow(data) - n.na + ess )
# 	
# 		# LL <- sum( lgamma(prior+counts) - lgamma(prior) )
# 		LL <- sum(lgamma(prior+counts)) - prod.sizes*lgamma(ess/prod.sizes) # faster!
# 		alpha_ij <- rowSums( prior, dims=length(pa) )
# 		N_ij <- rowSums( counts, dims=length(pa) )
# 		# return( LL + sum( lgamma(alpha_ij) - lgamma(alpha_ij+N_ij) ) )
# 		return( LL + prod.sizes.pa*lgamma(ess/prod.sizes.pa) - sum(lgamma(alpha_ij+N_ij)) ) # faster !
# 	}
# 	else
# 	{
# 		# correct for NAs with maximum a posteriori estimate
# 		n.na <- sum(is.na(data[,ni])) 
# 		counts <- counts + n.na * (counts + prior) / ( nrow(data) - n.na + ess )
# 	
# 		# LL <- sum( lgamma(prior+counts) - lgamma(prior) )
# 		LL <- sum(lgamma(prior+counts)) - prod.sizes*lgamma(ess/prod.sizes) # faster!
# 		return( LL - lgamma(1+sum(counts)) )
# 	}
# }

# compute.counts <- function( data, bins )
# {
# 	l <- length(bins)
# 	if( l > 1 ) # at least one parent
# 	{
# 		prev.cum.size <- c( 1, cumprod(bins[seq_len(l-1)]) )
# 		index <- data %*% prev.cum.size - sum(prev.cum.size) + 1
# 	}
# 	else
# 		index <- data
# 	counts <- tabulate( index, prod(bins) )
# 	dim(counts) <- bins
# 	return(counts)
# }

# find.best.parents <- function( aflml )
# {
# 	n.nodes <- nrow( aflml )
# 	bps <- matrix( 0, nrow=n.nodes, ncol=2^n.nodes )
# 	bss <- matrix( 0, nrow=n.nodes, ncol=2^n.nodes )
# 	bitmask <- 2^(seq_len(n.nodes)-1)
# 	
# 	for( ni in seq_len(n.nodes) )
# 		for( si in seq_len(2^n.nodes) ) # loop over all parent sets
# 		{
# 			if( bitAnd( si-1, bitmask[ni] ) ) # can't be parent of itself
# 				next
# 			
# 			bps[ni,si] <- si
# 			bss[ni,si] <- aflml[ni,si]
# 			
# 			# loop over all parents subsets that vary by 1 bit
# 			for( bi in which(bitAnd( si-1, bitmask ) > 0) )
# 			{
# 				ssi <- bitXor( si-1, bitmask[bi] ) + 1
# 				# cat(ni,si,bss[ni,si],bitmask,"\n")
# 				# cat(ni,ssi,bss[ni,ssi],bitmask,"\n")
# 				if( bss[ni,ssi] > bss[ni,si] )
# 				{
# 					bps[ni,si] <- bps[ni,ssi]
# 					bss[ni,si] <- bss[ni,ssi]
# 				}
# 			}	
# 		}
# 	return(bps)	
# }

# find.best.sinks <- function( bps, aflml )
# {
# 	n.nodes <- nrow(bps)
# 	bitmask <- 2^(seq_len(n.nodes)-1)
# 
# 	scores <- rep( 0, 2^n.nodes )
# 	sinks <- rep( -1, 2^n.nodes )
# 	
# 	for( si in seq_len(2^n.nodes) )
# 		for( sink in which(bitAnd( si-1, bitmask ) > 0) )
# 		{
# 			upvars <- bitXor( si-1, bitmask[sink] ) + 1
# 			skore <- scores[upvars] + aflml[sink, bps[sink, upvars]]
# 			
# 			if( sinks[si] == -1 || skore > scores[si] )
# 			{
# 				scores[si] <- skore
# 				sinks[si] <- sink
# 			}	
# 		}
# 	return(sinks)
# }

find.best.ordering <- function( sinks )
{
	n.nodes <- log2(length(sinks))
	bitmask <- 2^(seq_len(n.nodes)-1)

	order <- rep( 0, n.nodes )
	left <- 2^n.nodes

	for( i in n.nodes:1 )
	{
   	order[i] <- sinks[left]
   	left <- bitAnd( left-1, bitFlip(bitmask[order[i]], bitWidth=n.nodes) ) + 1
	}
	
	return( order )
}

find.best.network <- function( order, bps )
{
	n.nodes <- length(order)
	bitmask <- 2^(seq_len(n.nodes)-1)
	parents <- rep( 0, n.nodes )
	predecs <- 1

	for( i in seq_len(n.nodes) )
	{
		parents[order[i]] <- bps[order[i], predecs];
		predecs <- bitOr(predecs-1, bitmask[order[i]] ) + 1;
	}
	
	return( parents )
}

# convert array of parent sets to a full dag
parents.to.dag <- function( parents )
{
	n.nodes <- length( parents );
	bitmask <- 2^(seq_len(n.nodes)-1)
	dag <- matrix( 0, n.nodes, n.nodes);

	for( ni in seq_len(n.nodes) )
	{
   	bits <- which( bitAnd(parents[ni]-1,bitmask) > 0 );
    	dag[bits, ni] <- 1;
	}
	
	return( dag )
}

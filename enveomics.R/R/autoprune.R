
enve.prune.dist <- function
### Automatically prunes a tree, to keep representatives of each clade.
   (t,
### A `phylo` object or a path to the Newick file.
   dist.quantile=0.25,
### The quantile of edge lengths.
   min_dist,
### The minimum distance to allow between two tips. If not set, dist.quantile is
### used instead to calculate it.
   quiet=FALSE,
### Boolean indicating if the function must run without output.
   max_iters=100,
### Maximum number of iterations.
   min_nodes_random=4e4,
### Minimum number of nodes to trigger "tip-pairs" nodes sampling. This sampling
### is less reproducible and more computationally expensive, but it's the only
### solution if the cophenetic matrix exceeds 2^31-1 entries; above that, it
### cannot be represented in R.
   random_nodes_frx=1
### Fraction of the nodes to be sampled if more than `min_nodes_random`.
   ){
   if(!requireNamespace("ape", quietly=TRUE))
      stop('Unavailable ape library.');
   if(is.character(t)) t <- ape::read.tree(t)
   if(missing(min_dist)){
      if(dist.quantile>0){
	 min_dist <- as.numeric(quantile(t$edge.length, dist.quantile));
      }else{
         min_dist <- as.numeric(min(t$edge.length[t$edge.length>0]));
      }
   }
   if(!quiet) cat('\nObjective minimum distance: ',min_dist,', initial tips: ',length(t$tip.label),'\n', sep='');
   round=1;
   while(round <= max_iters){
      if(length(t$tip.label) > min_nodes_random){
	 if(!quiet) cat('  | Iter: ',round-1,', Tips: ', length(t$tip.label),
	 	', reducing tip-pairs.\n', sep='');
         rnd.nodes <- sample(t$tip.label, length(t$tip.label)*random_nodes_frx);
	 t <- enve.__prune.reduce(t, rnd.nodes, min_dist, quiet);
      }else{
	 if(!quiet) cat(' Gathering distances...\r');
	 d <- cophenetic(t);
	 diag(d) <- NA;
	 if(!quiet) cat('  | Iter: ',round-1,', Tips: ', length(t$tip.label),
		', Median distance: ', median(d, na.rm=TRUE),
      		', Minimum distance: ', min(d, na.rm=TRUE),
		'\n', sep='');
	 # Run iteration
	 if(min(d, na.rm=TRUE) < min_dist){
	    t <- enve.__prune.iter(t, d, min_dist, quiet);
	 }else{
	    break;
	 }
      }
      round <- round + 1;
   }
   return(t);
### Returns a pruned phylo object.
}

enve.__prune.reduce <- function
### Internal function for enve.prune.dist
   (t, nodes, min_dist, quiet){
   if(!quiet) pb <- txtProgressBar(1, length(nodes), style=3);
   for(i in 1:length(nodes)){
      node.name <- nodes[i];
      if(!quiet) setTxtProgressBar(pb, i);
      # Get node ID
      node <- which(t$tip.label==node.name);
      if(length(node)==0) next;
      # Get parent and distance to parent
      parent.node <- t$edge[ t$edge[,2]==node, 1];
      # Get edges to parent
      parent.edges <- which(t$edge[,1]==parent.node);
      stopit <- FALSE;
      for(j in parent.edges){
	 for(k in parent.edges){
	    if(j != k & t$edge[j,2]<length(t$tip.label) & t$edge[k,2]<length(t$tip.label) & sum(t$edge.length[c(j,k)]) < min_dist){
	       t <- ape::drop.tip(t, t$edge[k,2]);
	       stopit <- TRUE;
	       break;
	    }
	 }
	 if(stopit) break;
      }
   }
   if(!quiet) cat('\n');
   return(t);
}

enve.__prune.iter <- function
### Internal function for enve.prune.dist
   (t,
   dist,
   min_dist,
   quiet){
   ori_len <- length(t$tip.label);
   # Prune
   if(!quiet) pb <- txtProgressBar(1, ncol(dist)-1, style=3);
   ignore <- c();
   for(i in 1:(ncol(dist)-1)){
      if(i %in% ignore) next;
      for(j in (i+1):nrow(dist)){
	 if(dist[j, i]<min_dist){
	    t <- ape::drop.tip(t, rownames(dist)[j]);
	    ignore <- c(ignore, j);
	    break;
	 }
      }
      if(!quiet) setTxtProgressBar(pb, i);
   }
   if(!quiet) cat('\n');
   # Check if it droped tips
   cur_len <- length(t$tip.label);
   if(cur_len == ori_len){
      stop("Internal error: small edge found in tree, with no equivalent in distance matrix.\n");
   }
   return(t);
}


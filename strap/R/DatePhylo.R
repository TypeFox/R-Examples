DatePhylo <- function(tree, ages, rlen=0, method="basic", add.terminal=FALSE) {

  # Stop if using Ruta method but not supplying a tree with branch lengths:
  if(is.null(tree$edge.length) && method == "ruta") stop("Tree has no branch lengths (required for Ruta method).")

  # Stop if taxon names do not match between tree and ages matrix:
  if(length(c(setdiff(tree$tip.label, rownames(ages)), setdiff(rownames(ages), tree$tip.label))) > 0) stop("Taxon names for tree and ages do not match.")

  # Stop if the ages matrix does not have columns namd "FAD" and "LAD":
  if(length(match(c("FAD", "LAD"), colnames(ages))) < 2) stop("Ages matrix must have FAD and LAD (First and Last Appearance Datum) columns.")

  # Stop if any FAD is younger than its LAD:
  if(any(ages[, "FAD"] < ages[, "LAD"])) stop("FADs must all be at least as old as LADs.")

  # Stop if root length is negative:
  if(rlen < 0) stop("Root length cannot be negative.")

  # Stop if problem with root length:
  if(rlen == 0 && method != "basic") stop("If not using the basic method then rlen (root length) must be a positive value.")

  # Stop if method not available:
  if(method != "basic" && method != "ruta" && method != "equal") stop("Method must be one of basic, equal, or ruta.")

  # If method is not Ruta set all branch lengths equal to 1:
  if(is.null(tree$edge.length) || method != "ruta") tree$edge.length <- rep(1, length(tree$edge[, 1]))

  # Get node numbers to start:
  nodes <- c(Ntip(tree) + 1):(Nnode(tree) + Ntip(tree))
  
  # Create vector for node ages:
  node.ages <- nodes
  
  # Get starting node ages:
  for (i in 1:length(nodes)) node.ages[i] <- max(ages[tree$tip.label[FindDescendants(nodes[i], tree)], "FAD"])
  
  # Get all ages (tips and nodes):
  all.ages <- as.vector(c(ages[tree$tip.label, "FAD"], node.ages))
  
  # Set root age to maximum age plus root length:
  all.ages[Ntip(tree) + 1] <- all.ages[Ntip(tree) + 1] + rlen
  
  # Create time-scaled tree:
  time.tree <- tree
  
  # Add branch lengths as time:
  time.tree$edge.length <- abs(apply(matrix(all.ages[tree$edge], ncol=2), 1, diff))

  # Only continue if non basic dating option chosen:
  if(method != "basic") {
    
    # Keep going until there are no zero-length branches:
    while(min(time.tree$edge.length[grep(TRUE, tree$edge.length > 0)]) == 0) {
      
      # Record top zero-length branch encountered that is not also a zero change branch if using the Ruta method:
      share.branches <- intersect(grep(TRUE, time.tree$edge.length == 0), grep(TRUE, tree$edge.length > 0))[order(dist.nodes(tree)[(Ntip(tree) + 1), ][tree$edge[intersect(grep(TRUE, time.tree$edge.length == 0), grep(TRUE, tree$edge.length > 0)), 2]], decreasing=TRUE)[1]]
      
      # Keep going until there is a positive length branch:
      while(max(time.tree$edge.length[share.branches]) == 0) {
        
        # Find branches ancestral to those in memory:
        share.branches <- unique(c(share.branches, match(time.tree$edge[share.branches, 1], time.tree$edge[, 2])))

      }

      # Get total branch time:
      branch.time <- sum(time.tree$edge.length[share.branches])
      
      # Get number of branches to share:
      n.branches.to.share <- length(share.branches)
      
      # Get novel node ages (based on equal method):
      new.node.ages <- seq(from=range(all.ages[time.tree$edge[share.branches]])[1], to=range(all.ages[time.tree$edge[share.branches]])[2], by=sum(time.tree$edge.length[share.branches]) / n.branches.to.share)

      # Case if dating method is Ruta:
      if(method == "ruta") {
        
        # Get branch proportion based on branch lengths from input tree:
        branch.proportions <- ((tree$edge.length[share.branches] / sum(tree$edge.length[share.branches])))
        
        # Update novel node ages:
        new.node.ages <- c(new.node.ages[1], cumsum(branch.proportions[1:(length(branch.proportions) - 1)] * diff(range(new.node.ages))) + min(range(new.node.ages)), new.node.ages[length(new.node.ages)])

      }
      
      # Update node ages:
      all.ages[unique(as.vector(time.tree$edge[share.branches, 2:1]))] <- new.node.ages

      # Update branch lengths as time:
      time.tree$edge.length <- abs(apply(matrix(all.ages[tree$edge], ncol=2), 1, diff))

    }

  }
  
  # Add ranges of taxa to terminal branch lengths if requested:
  if(add.terminal) time.tree$edge.length[match(1:Ntip(time.tree), time.tree$edge[, 2])] <- time.tree$edge.length[match(1:Ntip(time.tree), time.tree$edge[, 2])] + ages[time.tree$tip.label, "FAD"] - ages[time.tree$tip.label, "LAD"]

  # Get root age:
  root.age <- max(all.ages)

  # Store root age:
  time.tree$root.time <- root.age
  
  # Output tree:
  return(time.tree)

}
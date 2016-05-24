StratPhyloCongruence <- function(trees, ages, rlen=0, method="basic", samp.perm=1000, rand.perm=1000, hard=TRUE, randomly.sample.ages=FALSE, fix.topology=TRUE, fix.outgroup=TRUE) {

  # Make sure permutation numbers are not negative:
  if(samp.perm <= 0 || rand.perm <= 0) stop("Number of permutations must be positive.")

  # Calculate sum of stratigraphic ranges of taxa (SRL):
  SRL <- sum(ages[, "FAD"] - ages[, "LAD"])
  
  # If sum of stratigraphic ranges is zero (all FADs are equal to their LADs):
  if(SRL == 0) {
  
    # If user has requested random sampling of ages method:
    if(randomly.sample.ages == TRUE) {
    
      # Warn user:
      cat("Can not perform random sampling of ages method as sum of stratigraphic ranges is zero.\nFor RCI sum of stratigraphic ranges is set to one.\n\n")
      
      # Set randomly sampled ages as false:
      randomly.sample.ages <- FALSE
      
    # If randomly sampled ages is set to false:
    } else {

      # Warn user:
      cat("Sum of stratigraphic ranges is set to one (to avoid dividing by zero).\n\n")
    
    }
    
    # Set SRL as one:
    SRL <- 1
    
  }
  
  # Create separate ages matrices for sampled and random permutations:
  samp.ages <- rand.ages <- ages

  # Set as null variables that may not be used later but are still included in the output:
  samp.permutations <- samp.trees <- NULL
  
  # Set root age:
  root.age <- max(ages[, "FAD"]) + rlen
  
  # If a single tree:
  if(class(trees) == "phylo") {
    
    # Convert to list to ensure a standard format for input tree(s):
    trees <- list(trees)
    
    # Set class as multiPhylo:
    class(trees) <- "multiPhylo"

  }

  # Get outgroup taxon for later re-rooting:
  if(fix.outgroup) {

    # Find maximum node height (indicating rooting taxon):
    max.node.height<- max(node.height(trees[[1]])[1:Ntip(trees[[1]])])

    # Find rows with maximum node height:
    max.node.height.rows <- grep(TRUE, node.height(trees[[1]])[1:Ntip(trees[[1]])] == max.node.height)
    
    # If there is not a single outgroup taxon for rooting:
    if(length(max.node.height.rows) > 1) {

      # Reset fix.outgroup to FALSE:
      fix.outgroup <-  FALSE

      # Warn user:
      print("Tree doesn't have single outgroup, fix.outgroup set to FALSE.")

    } else {

      # Set outgroup taxon:
      outgroup.taxon <- trees[[1]]$tip.label[max.node.height.rows]

    }
  
  }

  # Get total number of input trees:
  Ntrees <- length(trees)
  
  # If trees have edge lengths and are using Ruta dating method:
  if(length(grep("edge.length", names(unlist(trees)))) > 0 && method == "ruta") {
  
    # Collect all edge lengths in a single vector:
    edge.lengths <- as.numeric(unlist(trees)[grep("edge.length", names(unlist(trees)))])
  
  # If trees do not have edge lengths or if not using Ruta method:
  } else {
    
    # Set edge lengths vector to NULL:
    edge.lengths <- NULL

  }
  
  # If user wants to fix the topology as in GERt:
  if(fix.topology) {
    
    # Create random trees by resampling input trees:
    rand.trees <- trees[sample(1:length(trees), rand.perm, replace=TRUE)]
    
    # Reset class as multiPhylo:
    class(rand.trees) <- "multiPhylo"
    
    # For each random tree:
    for(i in 1:length(rand.trees)) {
      
      # Make sure it is fully bifurcating:
      rand.trees[[i]] <- multi2di(rand.trees[[i]])
      
      # If fixing the outgroup taxon (recommended):
      if(fix.outgroup) {

        # Randomly sample taxon names:
        sampled.names <- sample(rand.trees[[i]]$tip.label)

        # Ensure outgroup taxon is listed first:
        sampled.names <- c(outgroup.taxon, sampled.names[grep(FALSE, sampled.names == outgroup.taxon)])

        # Overwrite tip names with sampeld names:
        rand.trees[[i]]$tip.label <- sampled.names

      # If not fixing the outgroup taxon:
      } else {
  
        # Shuffle all tip labels:
        rand.trees[[i]]$tip.label <- sample(rand.trees[[i]]$tip.label)

      }
  
    }
    
  # If user does not want to fix the topology:
  } else {
  
    # Make random trees for permutations:
    rand.trees <- rmtree(rand.perm, Ntip(trees[[1]]))
  
  }

  # If edge.lengths is not NULL:
  if(!is.null(edge.lengths)) {
    
    # For each random tree:
    for(i in 1:rand.perm) {
      
      # Draw edge lengths from a uniform distribution between the limits of observed edge lengths:
      rand.trees[[i]]$edge.length <- runif(length(rand.trees[[i]]$edge.length), min=min(edge.lengths), max=max(edge.lengths))
    
    }

  }

  # For each tree in trees:
  for(i in 1:length(trees)) {
  
    # Date tree using settings specified in function:
    trees[[i]] <- DatePhylo(tree=trees[[i]], ages=ages, rlen=rlen, method=method, add.terminal=FALSE)
  
  }
  
  # If using the sample permutation loop:
  if(hard == FALSE || randomly.sample.ages == TRUE) {

    # Randomly sample numbers from 1 to Ntrees:
    tree.nos <- sample(c(1:Ntrees), samp.perm, replace=TRUE)
    
    # Create sampled trees list:
    samp.trees <- rmtree(samp.perm, Ntip(trees[[1]]))
    
    # For each sample permutation:
    for(i in 1:samp.perm) {
    
      # If polytomies are to be considered hard:
      if(hard == TRUE) {
        
        # Store randomly sampled tree:
        samp.trees[[i]] <- trees[[tree.nos[i]]]
      
      # If polytomies are to be considered soft:
      } else {
        
        # Store randomly sampled randomly bifurcated tree:
        samp.trees[[i]] <- multi2di(trees[[tree.nos[i]]])

      }
      
    }
    
  # If not permuting ages of polytomie resolutions:
  } else {
  
    # Reset samp.perm to zero:
    samp.perm <- 0

  }

  # Additional steps if randomly sampling ages:
  if(randomly.sample.ages) {
    
    # Create matrix to store randomly sampled ages ages:
    randomly.sampled.ages <- matrix(nrow=nrow(ages), ncol=max(c(samp.perm, rand.perm)) * 2)
    
    # For each taxon:
    for(i in 1:nrow(ages)) {
    
      # Draw random ages from a uniform distribution between FAD and LAD:
      randomly.sampled.ages[i, ] <- runif(samp.perm * 2, min=ages[i, "LAD"], ages[i, "FAD"])

    }
    
    # Take first half of matrix:
    firsthalf <- as.vector(randomly.sampled.ages[, 1:(length(randomly.sampled.ages[1,]) / 2)])
    
    # Take second half of matrix:
    secondhalf <- as.vector(randomly.sampled.ages[, ((length(randomly.sampled.ages[1,]) / 2) + 1):length(randomly.sampled.ages[1,])])
    
    # Sort halves relative to each other to get FADs and LADs:
    randomly.sampled.ages <- apply(cbind(firsthalf, secondhalf), 1, sort)
    
    # Define randomly sampled ages FADs:
    randomly.sampled.ages.FAD <- matrix(randomly.sampled.ages[2,], nrow=nrow(ages))
    
    # Define randomly sampled ages LADs:
    randomly.sampled.ages.LAD <- matrix(randomly.sampled.ages[1,], nrow=nrow(ages))

    # Add rownames:
    rownames(randomly.sampled.ages.FAD) <- rownames(randomly.sampled.ages.LAD) <- rownames(ages)
    
    # Set randomly sampled ages root lengths:
    randomly.sampled.ages.rlen <- rlen + max(ages) - apply(randomly.sampled.ages.FAD, 2, max)
    
    # Calculate Gmax (the worst possible fit to stratigraphy):
    randomly.sampled.ages.Gmax <- apply(root.age - randomly.sampled.ages.FAD, 2, sum)
    
    # Calculate Gmin (the best possible fit to stratigraphy - if no ADRs):
    randomly.sampled.ages.Gmin <- (apply(apply(apply(randomly.sampled.ages.FAD, 2, sort), 2, diff), 2, sum)) + (2 * (root.age - apply(randomly.sampled.ages.FAD, 2, max)))
    
    # Calculate sum of stratigraphic ranges of taxa (SRL):
    randomly.sampled.ages.SRL <- apply(randomly.sampled.ages.FAD - randomly.sampled.ages.LAD, 2, sum)
    
  }
  
  # Create matrix to store output with extra column for GERt:
  input.permutations <- matrix(ncol=10, nrow=length(trees))
  
  # Annotate matrix:
  colnames(input.permutations) <- c("SCI", "RCI", "GER", "MSM*", "est.p.SCI", "est.p.RCI", "est.p.GER", "est.p.MSM*", "GER*", "GERt")
  rownames(input.permutations) <- paste("tree_",c(1:length(trees)),sep="")  

  # Calculate Gmax (the worst possible fit to stratigraphy):
  Gmax <- sum(root.age - ages[, "FAD"])
  
  # Calculate Gmin (the best possible fit to stratigraphy - if no ADRs):
  Gmin <- diff(range(ages[, "FAD"])) + (2 * (root.age - max(ages[, "FAD"])))

  # Start progress bar:
  pb <- txtProgressBar(min = 0, max = (Ntrees + samp.perm + rand.perm), style = 3)

  # Text bar counter:
  counter <- 0
  
  # Vector to store random MIGs for calculating GER* and GERt:
  rand.MIG <- samp.MIG <- obs.MIG <- vector(mode="numeric")

  # Cycle through each tree in the list:
  for (i in 1:length(trees)) {
  
    # Isolate ith tree:
    tree <- trees[[i]]
    
    # MIG calculation:
    obs.MIG[i] <- MIG <- sum(tree$edge.length)
    
    # Create node matrix for SCI calculation:
    newnodes <- matrix(ncol=2, nrow=Nnode(tree) - 1)
    
    # Add rownames to matrix:
    rownames(newnodes) <- c((Ntip(tree) + 2):(Ntip(tree) + Nnode(tree)))
    
    # Add column names to matrix:
    colnames(newnodes) <- c("anc", "dec")
    
    # For each node (excluding the root):
    for (j in 1:nrow(newnodes)) {
      
      # Get descendants for current node:
      dec <- tree$tip.label[FindDescendants(as.numeric(rownames(newnodes)[j]), tree)]
      
      # Get descendants for immediate ancestral node:
      anc <- tree$tip.label[FindDescendants(tree$edge[match(as.numeric(rownames(newnodes)[j]), tree$edge[, 2]), 1], tree)]
      
      # Remove descendnats from ancestral node already found in current node:
      anc <- anc[match(anc, dec, nomatch=0) == 0]
      
      # Store oldest age for current node descendants:
      newnodes[j, 2] <- max(ages[dec, "FAD"])
      
      # Store oldest age for ancestral node descendants:
      newnodes[j, 1] <- max(ages[anc, "FAD"])
      
    }
    
    # SCI calculation:
    input.permutations[i, "SCI"] <- sum(newnodes[, "anc"] >= newnodes[, "dec"]) / nrow(newnodes)
    
    # RCI calculation:
    input.permutations[i, "RCI"] <- (1 - (MIG / SRL)) * 100
    
    # GER calculation:
    input.permutations[i, "GER"] <- 1 - ((MIG - Gmin) / (Gmax - Gmin))
    
    # MSM* calculation:
    input.permutations[i, "MSM*"] <- Gmin / MIG
    
    # Update counter:
    counter <- counter + 1
    
    # Update progress bar:
    setTxtProgressBar(pb, counter)
    
  }
  
  # If polytomies are not hard (i.e. if they are soft):
  if(hard == FALSE || randomly.sample.ages == TRUE) {
    
    # Create matrix to store permutation results with extra column for GERt:
    samp.permutations <- matrix(nrow=samp.perm, ncol=10)

    # Add column names:
    colnames(samp.permutations) <- c("SCI", "RCI", "GER", "MSM*", "est.p.SCI", "est.p.RCI", "est.p.GER", "est.p.MSM*", "GER*", "GERt")

    # For each sample permutation:
    for(i in 1:samp.perm) {
      
      # Isolate ith sample tree:
      samp.tree <- samp.trees[[i]]
      
      # If using randomly sampled ages:
      if(randomly.sample.ages) {
        
        # Overwrite rand.ages with ith randomly.sampled.ages
        samp.ages[, "FAD"] <- randomly.sampled.ages.FAD[rownames(rand.ages), i]
        
        # Overwrite rand.ages with ith randomly.sampled.ages
        samp.ages[, "LAD"] <- randomly.sampled.ages.LAD[rownames(rand.ages), i]
        
        # Overwrite rlen with ith randomly sampled ages rlen:
        rlen <- randomly.sampled.ages.rlen[i]
        
        # Overwrite Gmin with randomly sampled ages Gmin:
        Gmin <- randomly.sampled.ages.Gmin[i]
        
        # Overwrite Gmax with randomly sampled ages Gmax:
        Gmax <- randomly.sampled.ages.Gmax[i]
        
        # Overwrite SRL with randomly sampled ages SRL:
        SRL <- randomly.sampled.ages.SRL[i]
        
      }
      
      # Date tree using settings specified in function:
      samp.tree <- DatePhylo(tree=samp.tree, ages=samp.ages, rlen=rlen, method=method, add.terminal=FALSE)
      
      # MIG calculation:
      samp.MIG[i] <- MIG <- sum(samp.tree$edge.length)
      
      # Create node matrix for SCI calculation:
      newnodes <- matrix(ncol=2, nrow=Nnode(samp.tree) - 1)
      
      # Add rownames to matrix:
      rownames(newnodes) <- c((Ntip(samp.tree) + 2):(Ntip(samp.tree) + Nnode(samp.tree)))
      
      # Add column names to matrix:
      colnames(newnodes) <- c("anc", "dec")
      
      # For each node (excluding the root):
      for (j in 1:nrow(newnodes)) {
        
        # Get descendants for current node:
        dec <- samp.tree$tip.label[FindDescendants(as.numeric(rownames(newnodes)[j]), samp.tree)]
        
        # Get descendants for immediate ancestral node:
        anc <- samp.tree$tip.label[FindDescendants(samp.tree$edge[match(as.numeric(rownames(newnodes)[j]), samp.tree$edge[, 2]), 1], samp.tree)]
        
        # Remove descendnats from ancestral node already found in current node:
        anc <- anc[match(anc, dec, nomatch=0) == 0]
        
        # Store oldest age for current node descendants:
        newnodes[j, 2] <- max(samp.ages[dec, "FAD"])
        
        # Store oldest age for ancestral node descendants:
        newnodes[j, 1] <- max(samp.ages[anc, "FAD"])
        
      }
      
      # SCI calculation:
      samp.permutations[i, "SCI"] <- sum(newnodes[, "anc"] >= newnodes[, "dec"]) / nrow(newnodes)
      
      # RCI calculation:
      samp.permutations[i, "RCI"] <- (1 - (MIG / SRL)) * 100
      
      # GER calculation:
      samp.permutations[i, "GER"] <- 1 - ((MIG - Gmin) / (Gmax - Gmin))
      
      # MSM* calculation:
      samp.permutations[i, "MSM*"] <- Gmin / MIG

      # Save sampled tree:
      samp.trees[[i]] <- samp.tree
      
      # Update counter:
      counter <- counter + 1
      
      # Update progress bar:
      setTxtProgressBar(pb, counter)
      
    }

    # Add MIG to the sample permutaion results:
    samp.permutations <- cbind(samp.permutations, samp.MIG)

    # Add column heading for MIG:
    colnames(samp.permutations)[length(colnames(samp.permutations))] <- "MIG"

  }

  # Create matrix to store permutation results:
  rand.permutations <- matrix(nrow=rand.perm, ncol=4)
  
  # Add column names:
  colnames(rand.permutations) <- c("SCI", "RCI", "GER", "MSM*")
  
  # For each permutation:
  for (i in 1:rand.perm) {
    
    # Isolate ith random tree:
    rand.tree <- rand.trees[[i]]

    # If fixing the outgroup (recommended):
    if(fix.outgroup == TRUE && fix.topology == FALSE) {

      # Overwrite random names with real taxon names
      rand.tree$tip.label <- sample(trees[[1]]$tip.label)

      # Re-root tree on outgroup taxon:
      rand.tree <- root(rand.tree, outgroup.taxon, resolve.root=TRUE)
    
    }

    # If not fixing the outgroup:
    if(fix.outgroup == FALSE) {

      # Overwrite random names with real taxon names
      rand.tree$tip.label <- sample(trees[[1]]$tip.label)

    }

    # If using randomly sampled ages:
    if(randomly.sample.ages) {
      
      # Overwrite rand.ages with ith randomly.sampled.ages
      rand.ages[, "FAD"] <- randomly.sampled.ages.FAD[rownames(rand.ages), i]
      
      # Overwrite rand.ages with ith randomly.sampled.ages
      rand.ages[, "LAD"] <- randomly.sampled.ages.LAD[rownames(rand.ages), i]

      # Overwrite rlen with ith randomly sampled ages rlen:
      rlen <- randomly.sampled.ages.rlen[i]
      
      # Overwrite Gmin with randomly sampled ages Gmin:
      Gmin <- randomly.sampled.ages.Gmin[i]

      # Overwrite Gmax with randomly sampled ages Gmax:
      Gmax <- randomly.sampled.ages.Gmax[i]
      
      # Overwrite SRL with randomly sampled ages SRL:
      SRL <- randomly.sampled.ages.SRL[i]
      
    }
    
    # Date tree using settings specified in function:
    rand.tree <- DatePhylo(tree=rand.tree, ages=rand.ages, rlen=rlen, method=method, add.terminal=FALSE)
    
    # MIG calculation:
    rand.MIG[i] <- MIG <- sum(rand.tree$edge.length)
    
    # Create node matrix for SCI calculation:
    newnodes <- matrix(ncol=2, nrow=Nnode(rand.tree) - 1)
    
    # Add rownames to matrix:
    rownames(newnodes) <- c((Ntip(rand.tree) + 2):(Ntip(rand.tree) + Nnode(rand.tree)))
    
    # Add column names to matrix:
    colnames(newnodes) <- c("anc", "dec")
    
    # For each node (excluding the root):
    for (j in 1:nrow(newnodes)) {
      
      # Get descendants for current node:
      dec <- rand.tree$tip.label[FindDescendants(as.numeric(rownames(newnodes)[j]), rand.tree)]
      
      # Get descendants for immediate ancestral node:
      anc <- rand.tree$tip.label[FindDescendants(rand.tree$edge[match(as.numeric(rownames(newnodes)[j]), rand.tree$edge[, 2]), 1], rand.tree)]
      
      # Remove descendnats from ancestral node already found in current node:
      anc <- anc[match(anc, dec, nomatch=0) == 0]
      
      # Store oldest age for current node descendants:
      newnodes[j, 2] <- max(rand.ages[dec, "FAD"])
      
      # Store oldest age for ancestral node descendants:
      newnodes[j, 1] <- max(rand.ages[anc, "FAD"])
      
    }
    
    # SCI calculation:
    rand.permutations[i, "SCI"] <- sum(newnodes[, "anc"] >= newnodes[, "dec"]) / nrow(newnodes)
    
    # RCI calculation:
    rand.permutations[i, "RCI"] <- (1 - (MIG / SRL)) * 100
    
    # GER calculation:
    rand.permutations[i, "GER"] <- 1 - ((MIG - Gmin) / (Gmax - Gmin))
    
    # MSM* calculation:
    rand.permutations[i, "MSM*"] <- Gmin / MIG
    
    # Save randomisation tree:
    rand.trees[[i]] <- rand.tree
    
    # Update counter:
    counter <- counter + 1
    
    # Update progress bar:
    setTxtProgressBar(pb, counter)
  
  }
  
  # Close progress bar:
  close(pb)
  
  # Calculate input tree p.values for SCI:
  input.permutations[, "est.p.SCI"] <- pnorm(asin(sqrt(input.permutations[, "SCI"])), mean(asin(sqrt(rand.permutations[, "SCI"]))), sd(asin(sqrt(rand.permutations[, "SCI"]))), lower.tail=FALSE)

  # Calculate input tree p.values for RCI:
  input.permutations[, "est.p.RCI"] <- pnorm(input.permutations[, "RCI"], mean(rand.permutations[, "RCI"]), sd(rand.permutations[, "RCI"]), lower.tail=FALSE)
  
  # Calculate input tree p.values for GER:
  input.permutations[, "est.p.GER"] <- pnorm(asin(sqrt(input.permutations[, "GER"])), mean(asin(sqrt(rand.permutations[, "GER"]))), sd(asin(sqrt(rand.permutations[, "GER"]))), lower.tail=FALSE)
  
  # Calculate input tree p.values for MSM*:
  input.permutations[, "est.p.MSM*"] <- pnorm(asin(sqrt(input.permutations[, "MSM*"])), mean(asin(sqrt(rand.permutations[, "MSM*"]))), sd(asin(sqrt(rand.permutations[, "MSM*"]))), lower.tail=FALSE)
  
  # If sample permutations were performed:
  if(samp.perm > 0) {

    # Calculate sample permutation p.values for SCI:
    samp.permutations[, "est.p.SCI"] <- pnorm(asin(sqrt(samp.permutations[, "SCI"])), mean(asin(sqrt(rand.permutations[, "SCI"]))), sd(asin(sqrt(rand.permutations[, "SCI"]))), lower.tail=FALSE)
  
    # Calculate sample permutation p.values for RCI:
    samp.permutations[, "est.p.RCI"] <- pnorm(samp.permutations[, "RCI"], mean(rand.permutations[, "RCI"]), sd(rand.permutations[, "RCI"]), lower.tail=FALSE)
  
    # Calculate sample permutation p.values for GER:
    samp.permutations[, "est.p.GER"] <- pnorm(asin(sqrt(samp.permutations[, "GER"])), mean(asin(sqrt(rand.permutations[, "GER"]))), sd(asin(sqrt(rand.permutations[, "GER"]))), lower.tail=FALSE)
  
    # Calculate sample permutation p.values for MSM*:
    samp.permutations[, "est.p.MSM*"] <- pnorm(asin(sqrt(samp.permutations[, "MSM*"])), mean(asin(sqrt(rand.permutations[, "MSM*"]))), sd(asin(sqrt(rand.permutations[, "MSM*"]))), lower.tail=FALSE)
  
    # For each sampled tree:
    for(i in 1:samp.perm) {

      # Calculate GER* for ith sampled tree:
      samp.permutations[i, "GER*"] <- length(grep(TRUE, samp.MIG[i] <= rand.MIG)) / rand.perm

      # Calculate GERt:
      samp.permutations[i, "GERt"] <- 1 - ((samp.MIG[i] - min(rand.MIG)) / (max(rand.MIG) - min(rand.MIG)))
  
      # If GERt is negative re-scale to 0:
      if(samp.permutations[i, "GERt"] < 0) samp.permutations[i, "GERt"] <- 0
  
      # If GERt is greater than 1 re-scale to 1:
      if(samp.permutations[i, "GERt"] > 1) samp.permutations[i, "GERt"] <- 1

    }

  }

  # For each input tree:
  for (i in 1:Ntrees) {

    # Calculate GER* for ith input tree:
    input.permutations[i, "GER*"] <- length(grep(TRUE, obs.MIG[i] <= rand.MIG)) / rand.perm

    # Calculate GERt:
    input.permutations[i, "GERt"] <- 1 - ((obs.MIG[i] - min(rand.MIG)) / (max(rand.MIG) - min(rand.MIG)))

    # If GERt is negative re-scale to 0:
    if(input.permutations[i, "GERt"] < 0) input.permutations[i, "GERt"] <- 0

    # If GERt is greater than 1 re-scale to 1:
    if(input.permutations[i, "GERt"] > 1) input.permutations[i, "GERt"] <- 1
  
  }

  # Add MIG to input permutations:
  input.permutations <- cbind(input.permutations, obs.MIG)

  # Add MIG to random permutations:
  rand.permutations <- cbind(rand.permutations, rand.MIG)

  # Add MIG to column headings:
  colnames(input.permutations)[length(colnames(input.permutations))] <- "MIG"

  # Add MIG to column headings:
  colnames(rand.permutations)[length(colnames(rand.permutations))] <- "MIG"

  # If sample permutations were performed:
  if(samp.perm > 0) {

    # Add Wills p column:
    samp.permutations <- cbind(samp.permutations, rep(0, samp.perm))

    # Add Wills p to column headings:
    colnames(samp.permutations)[length(colnames(samp.permutations))] <- "p.Wills"

    # Add Wills p values:
    for(i in 1:samp.perm) samp.permutations[i, "p.Wills"] <- 1 - (sum(samp.permutations[i, "MIG"] < rand.permutations[, "MIG"]) / rand.perm)

  }

  # Add Wills p column:
  input.permutations <- cbind(input.permutations, rep(0, nrow(input.permutations)))

  # Add Wills p to column headings:
  colnames(input.permutations)[length(colnames(input.permutations))] <- "p.Wills"

  # Add Wills p values:
  for(i in 1:nrow(input.permutations)) input.permutations[i, "p.Wills"] <- 1 - (sum(input.permutations[i, "MIG"] < rand.permutations[, "MIG"]) / rand.perm)

  # Compile output as list:
  output <- list(input.permutations, samp.permutations, rand.permutations, trees, samp.trees, rand.trees)

  # Add names to output:
  names(output) <- c("input.tree.results", "samp.permutation.results", "rand.permutations", "input.trees", "samp.trees", "rand.trees")

  # Notify user that run is complete:
  cat(paste("Fit to stratigraphy measures for", nrow(input.permutations), "input trees have been calculated."))

  # Return output:
  return(invisible(output))


}
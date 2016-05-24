nested.trees <- function(xtree, xdata, Y.name, X.names, MaxTreeSize = NULL, family = "binomial", verbose = TRUE)
{
  if(!inherits(xtree, 'rpart')) stop('xtree have to be an rpart object!')
  ##	Nodes of the tree
  nodes <- as.numeric(rownames(xtree$frame))
  
  ##	Number of nodes and the last node
  nber_nodes <- length(nodes)
  last_node <- max(nodes)
  
  ##	Number of leaves
  nber_leaves <- sum(xtree$frame$var == "<leaf>")
  if(is.null(MaxTreeSize)) MaxTreeSize = nber_leaves
  
  ##	Leaves of the max-tree
  leaves_max_tree <- nodes[xtree$frame$var == "<leaf>"]
  
  nested_trees <- list()
  leaves <- list()
  deviances <- double()
  
  ##	Print 
  if(verbose) cat("Number of leaves in the max tree = ", nber_leaves, "\n")
  
  ##	Nested-tree with 2 leaves
  j = 1
  best_leaves <- leaves[[j]] <- c(2, 3)
  nested_trees[[j]] <- snip.rpart(xtree, toss = leaves[[j]])
  
  ##	Null model with only covariates in X
  fit_null <- glm(as.formula(paste(Y.name, " ~ ", paste(X.names, collapse = "+"))), data = xdata, family = family)
  
  deviance_null <- deviance(fit_null)
  
  fit_pltr_glm_given_tree <- tree2glm(nested_trees[[j]], xdata, Y.name, X.names, family = family)
  
  deviances[j] <- deviance(fit_pltr_glm_given_tree)
  
  
  nber <- 1 
  ##	Other nested-trees
  while(! setequal(best_leaves, leaves_max_tree))
  {
    nber <- nber + 1
    
    ##	Print
    if(verbose){
        cat("Best sub-tree \n", best_leaves, "\n")
        cat("\n Number of sub-trees = ", nber, "\n")
    }
    
    if(length(best_leaves) == MaxTreeSize)
    {
      cat("Max tree size ", MaxTreeSize, "reached \n")
      break
    }
    
    i <- 1
    List_leaves <- list()
    Vect_deviances <- c()
    for(l in best_leaves)
    {
      current_leaves <- c(setdiff(best_leaves, l), 2*l, 2*l+1)
      if(all(is.element(current_leaves, nodes)))
      {
        List_leaves[[i]] <- current_leaves
        current_tree = snip.rpart(xtree, toss = current_leaves)
        current_fit_pltr_lm_given_tree <- tree2glm(current_tree, xdata, Y.name, X.names, family = family)
        Vect_deviances[i] <- deviance(current_fit_pltr_lm_given_tree)
        i = i+1
      }
    }
    
    best_index <- which.min(Vect_deviances)
    best_current_leaves <- List_leaves[[best_index]]
    j = j + 1
    leaves[[j]] <- List_leaves[[best_index]]
    best_leaves <- List_leaves[[best_index]]
    deviances[j] <- Vect_deviances[best_index]
  }
  
  diff_deviances <- deviance_null - deviances
  
  return(list(leaves = leaves, null_deviance = deviance_null, deviances = deviances, diff_deviances = diff_deviances))
}

best.tree.BIC.AIC <- function(xtree, xdata, Y.name, X.names, family = "binomial", verbose = TRUE)
{
  time1 <- Sys.time()
  if(!inherits(xtree, 'rpart')) stop('xtree have to be an rpart object!')
  tree_size <- sum(xtree$frame$var == "<leaf>")    
  MaxTreeSize <- min(tree_size,10)
  
  nested_trees <- nested.trees(xtree, xdata, Y.name, X.names, MaxTreeSize = MaxTreeSize, family = family, verbose = verbose)
  dev <- c(nested_trees$null_deviance, nested_trees$deviances)
  
  # Compute BIC and AIC
  n <- nrow(xdata)
  
  BIC_AIC <- t(sapply(nested_trees$leaves, function(uu)
  {
    xxtree <- snip.rpart(xtree, toss = uu)
    
    xindicators_tree <- sapply(tree2indicators(xxtree), function(u) return(paste("as.integer(", u, ")")))
    
    xnber_indicators <- length(xindicators_tree)
    
    xformula <- as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "), "+", paste(xindicators_tree[-xnber_indicators], collapse = "+")))
    
    xupdate_lm <- glm(xformula, data = xdata, family = family)
    
    nber_parameters <- length(X.names) + length(uu)
    
    LL <- -deviance(xupdate_lm)
    
    return(c(LL - nber_parameters*log(n), LL - 2*nber_parameters))	
  }))
  
  ##	Add model with only X
  fit_null <- glm(as.formula(paste(Y.name, "~", paste(X.names, collapse = "+"))), data = xdata, family = family)
  
  nber_parameters <- length(X.names) + 1
  LL_null <- -deviance(fit_null)
  
  BIC_AIC <- rbind(c(LL_null - nber_parameters*log(n), LL_null - 2*nber_parameters), BIC_AIC)
  
  colnames(BIC_AIC) <- c("BIC", "AIC")
  
  best_index_BIC <- which.max(BIC_AIC[, "BIC"])
  best_index_AIC <- which.max(BIC_AIC[, "AIC"])
  
  nested_trees_leaves <- append(nested_trees$leaves, c(1), after = 0)
  
  best_tree_BIC <- snip.rpart(xtree, toss = nested_trees_leaves[[best_index_BIC]])
  best_tree_AIC <- snip.rpart(xtree, toss = nested_trees_leaves[[best_index_AIC]])
  
  ## Fit linear model for BIC
  if(best_index_BIC > 1)
  {
    indicators_tree <- sapply(tree2indicators(best_tree_BIC), function(u) return(paste("as.integer(", u, ")")))
    
    nber_indicators <- length(indicators_tree)
    
    xformula <- as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "), "+", paste(indicators_tree[-nber_indicators], collapse = "+")))
    
    lm_BIC <- glm(xformula, data = xdata, family = family)
  }else{
    lm_BIC <- glm(as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "))), data = xdata, family = family)
  }
  
  ## Fit linear model for AIC
  if(best_index_AIC > 1)
  {
    indicators_tree <- sapply(tree2indicators(best_tree_AIC), function(u) return(paste("as.integer(", u, ")")))
    
    nber_indicators <- length(indicators_tree)
    
    xformula <- as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "), "+", paste(indicators_tree[-nber_indicators], collapse = "+")))
    
    lm_AIC <- glm(xformula, data = xdata, family = family)
  }else{
    lm_AIC <- glm(as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "))), data = xdata, family = family)
  }
  ## Fit linear model for CV
  
  
  best_index <- c(BIC = best_index_BIC, AIC = best_index_AIC)
  
  best_tree <- list(BIC = best_tree_BIC, AIC = best_tree_AIC)
  
  fitLM <- list(BIC = lm_BIC, AIC = lm_AIC)
  
  time2 <- Sys.time()
  Timediff <- difftime(time2, time1)
  return(list(best_index = best_index, tree = best_tree, fit_glm = fitLM, Timediff = Timediff))
}

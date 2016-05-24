best.tree.permute <- function(xtree, xdata, Y.name, X.names, G.names, B = 10, args.rpart = list(cp = 0, minbucket = 20, maxdepth = 10), epsi = 1e-3, iterMax = 5, iterMin = 3, family = "binomial", LEVEL = 0.05, LB = FALSE, args.parallel = list(numWorkers = 1, type = "PSOCK"), verbose = TRUE)
{
  time1 <- Sys.time()
  if(!inherits(xtree, 'rpart')) stop('xtree have to be an rpart object!')
  ##	Fit null model
  fit_null = glm(as.formula(paste(Y.name, " ~ ", paste(X.names, collapse = "+"))), data = xdata, family = family)
  
  ##	Parameter of X covariates
  hat_gamma0 = fit_null$coef
  
  n = nrow(xdata)
  
  ##	Workers job
  MaxTreeList = list()
  Ind_0 <- which(xdata[X.names] == 0)
  Ind_1 <- which(xdata[X.names] == 1)
  Yp_list <- lapply(1:B, function(j) 
  {
    resp <- double(nrow(xdata))
    IndY_0 <- sample(Ind_0,length(Ind_0),replace = FALSE)
    IndY_1 <- sample(Ind_1,length(Ind_1),replace = FALSE)
    resp[Ind_0] <- xdata[Y.name][IndY_0,]
    resp[Ind_1] <- xdata[Y.name][IndY_1,]
    
    return(resp)
  })
  
  
  ##	Wrapper function to be executed in each machine
  wrapper <- function(xdata_p)
  {
    pltr_lm_p = pltr.glm(xdata_p, Y.name = "Yp", X.names = X.names, G.names = G.names, args.rpart = args.rpart, epsi = epsi, iterMax = iterMax, iterMin = iterMin, family = family, verbose = verbose)
    
    return(pltr_lm_p$tree)
  }
  
  ##	Parallel job
  numWorkers = args.parallel$numWorkers
  
  cat("\n ncores = ", numWorkers, " for bootstrap permutation !\n")
  
  List_xdatas <- lapply(Yp_list, function(Ypp) return(data.frame(Yp = Ypp, xdata)))
  
  
  MaxTreeList <- mclapply(List_xdatas, wrapper, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  
  ##	Max number of leaves "m"
  size_vect <- sapply(MaxTreeList, function(xxtree) {return(sum(xxtree$frame$var == "<leaf>"))})
  if (min(size_vect) <= 3)
  {
    supp_size <- which(size_vect <= 3)
    MaxTreeList[supp_size] <- NULL
    size_vect <- size_vect[-supp_size]
    Yp_list[supp_size] <- NULL
  }
  size_original_tree <- sum(xtree$frame$var == "<leaf>")
  
  m <- min(min(size_vect), size_original_tree, 16)
  B_adj <- length(size_vect)
  
  
  
  ##	Nested trees
  wrapper2 <- function(list_xtree_xdata)
  {
    return(nested.trees(xtree = list_xtree_xdata[[1]], xdata = list_xtree_xdata[[2]], Y.name = "Yp", X.names = X.names, MaxTreeSize = m, family = family, verbose = verbose)$diff_deviances)
  }
  
  List_xTrees_xDatas <- list()
  List_xTrees_xDatas <-	lapply(1:B_adj, function(j)
  {
    return(list(MaxTreeList[[j]], data.frame(Yp = Yp_list[[j]], xdata)))
  })
  
  
  cat("\n ncores = ", numWorkers, " for nested trees !\n")
  
  
  List_Diff_deviances <- mclapply(List_xTrees_xDatas, wrapper2, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  
  ##	Stop parallelization
  # 	stopCluster(cl)
  
  
  ##	Observed statistics
  obs_nested_trees = nested.trees(xtree, xdata, Y.name, X.names, MaxTreeSize = m, family = family, verbose = verbose)
  obs_diff_deviances = obs_nested_trees$diff_deviances
  
  ##	Compute bootstrap p.values
  Diff_deviances <- matrix(unlist(List_Diff_deviances), nrow = m - 1, byrow = FALSE)
  Diff_deviances <- cbind(obs_diff_deviances,Diff_deviances)
  mat_p_values <- t(apply(Diff_deviances,1,function(dev){sapply(dev,function(u){sum(dev >= u)/length(dev)})}))
  vect_p_values <- mat_p_values[,1]
  MinP <- apply(mat_p_values,2,min)
  p.val_selected <- sum(MinP <= MinP[1])/length(MinP)
  
  ##	Final selection
  index_best_tree <- which.min(vect_p_values)
  
  ##	Selected leaves
  selected_leaves <- obs_nested_trees$leaves[[index_best_tree]]
  
  ##	Selected tree
  selected_tree <- snip.rpart(xtree, toss = selected_leaves)
  
  ## Fit linear model for the selected tree
  
  indicators_tree = sapply(tree2indicators(selected_tree), function(u) return(paste("as.integer(", u, ")")))
  
  nber_indicators = length(indicators_tree)
  
  xformula = as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "), "+", paste(indicators_tree[-nber_indicators], collapse = "+")))
  
  selected_fit_glm = glm(xformula, data = xdata, family = family)
  
  selected_p_value <- vect_p_values[index_best_tree]
  
  Tree_Selected <- (p.val_selected < LEVEL)
  
  selected_fit_glm_tree <- list(fit_glm = selected_fit_glm, tree = selected_tree, p.value = selected_p_value, Tree_Selected = Tree_Selected)
  
  if(! Tree_Selected) 
  {
    selected_fit_glm = glm(as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "))), data = xdata, family = family)
  }
  
  time2 <- Sys.time()
  Timediff <- difftime(time2, time1)
  return(list(p.val_selected = p.val_selected, selected_model = selected_fit_glm_tree, fit_glm = selected_fit_glm, Timediff = Timediff, comp_p_values = vect_p_values, Badj = B_adj))
}

best.tree.bootstrap <- function(xtree, xdata, Y.name, X.names, G.names, B = 10, BB = 10, args.rpart = list(cp = 0, minbucket = 20, maxdepth = 10), epsi = 1e-3, iterMax = 5, iterMin = 3, family = "binomial", LEVEL = 0.05, LB = FALSE, args.parallel = list(numWorkers = 1), verbose = TRUE)
{
  time1<-Sys.time()
  if(!inherits(xtree, 'rpart')) stop('xtree have to be an rpart object!')
  ##  Fit null model
  fit_null <- glm(as.formula(paste(Y.name, " ~ ", paste(X.names, collapse = "+"))), data = xdata, family = family)
  
  ##  Null deviance
  #   null_deviance = deviance(fit_null)
  
  ##	Parameter of X covariates
  hat_gamma0 <- fit_null$coef
  
  n <- nrow(xdata)
  
  ##	For product
  product <- ifelse(length(X.names) == 1, "*", "%*%")
  
  ##	Workers job
  MaxTreeList <- list()
  
  
  if(family == "binomial")
  {
    Vect_PY1 <- fit_null$fitted
    Yb_list <- list()
    Yb_list <- lapply(1:B, function(j) 
    {
      return(sapply(Vect_PY1, function(xp) {
        return(sample(c(1, 0), 1, replace = TRUE, prob = c(xp, 1-xp)))
      }))
    })
  }else{
    df <- length(hat_gamma0)
    hat_sd2e <- sqrt(sum(fit_null$residuals^2)/(n-df))
    
    Yb_list <- list()
    Yb_list <- lapply(1:B, function(j)
    {
      return(fit_null$fitted.values + rnorm(n, mean = 0, sd = hat_sd2e))
    })
  }
  
  ##	Wrapper function to be executed in each machine
  wrapper <- function(xdata_b)
  {
    pltr_lm_b <- pltr.glm(xdata_b, Y.name = "Yb", X.names = X.names, G.names = G.names, args.rpart = args.rpart, epsi = epsi, iterMax = iterMax, iterMin = iterMin, family = family, verbose = verbose)
    
    return(pltr_lm_b$tree)
  }
  
  ##	Parallel job
  numWorkers <- args.parallel$numWorkers
  cat("\n ncores = ", numWorkers, " for  inner layer bootstrap !\n")
  
  List_xdatas <- lapply(Yb_list, function(Ybb) return(data.frame(Yb = Ybb, xdata)))
  
  MaxTreeList <- mclapply(List_xdatas, wrapper, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  
  
  size_vect <- sapply(MaxTreeList, function(xxtree) {return(sum(xxtree$frame$var == "<leaf>"))})
  if (min(size_vect) <= 3)
  {
    supp_size <- which(size_vect <= 3)
    MaxTreeList[supp_size] <- NULL
    size_vect <- size_vect[-supp_size]
    Yb_list[supp_size] <- NULL
  }
  size_original_tree <- sum(xtree$frame$var == "<leaf>")
  
  m <- min(min(size_vect), size_original_tree, 16)
  B_adj <- length(size_vect)  
  
  ##	Nested trees
  wrapper2 <- function(list_xtree_xdata)
  {
    return(nested.trees(xtree = list_xtree_xdata[[1]], xdata = list_xtree_xdata[[2]], Y.name = "Yb", X.names = X.names, MaxTreeSize = m, family = family, verbose = verbose)$diff_deviances)
  }
  
  List_xTrees_xDatas <- list()
  List_xTrees_xDatas <-	lapply(1:B_adj, function(j)
  {
    return(list(MaxTreeList[[j]], data.frame(Yb = Yb_list[[j]], xdata)))
  })
  
  
  cat("\n ncores = ", numWorkers, " for nested trees in the inner layer !\n")
  
  List_Diff_deviances <- mclapply(List_xTrees_xDatas, wrapper2, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  
  ##	Observed statistics
  obs_nested_trees <- nested.trees(xtree, xdata, Y.name, X.names, MaxTreeSize = m, family = family, verbose = verbose)
  obs_diff_deviances <- obs_nested_trees$diff_deviances
  
  ##	Compute bootstrap p.values
  Diff_deviances <- matrix(unlist(List_Diff_deviances), nrow = m - 1, byrow = FALSE)
  vect_p_values <- sapply(1:(m-1), function(j) return(mean(Diff_deviances[j,] > obs_diff_deviances[j])))
  
  ##	Final selection
  index_best_tree <- which.min(vect_p_values)
  
  ##	Selected leaves
  selected_leaves <- obs_nested_trees$leaves[[index_best_tree]]
  
  ##	Selected tree
  selected_tree <- snip.rpart(xtree, toss = selected_leaves)
  
  ## Fit linear model for the selected tree
  
  indicators_tree <- sapply(tree2indicators(selected_tree), function(u) return(paste("as.integer(", u, ")")))
  
  nber_indicators <- length(indicators_tree)
  
  xformula <- as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "), "+", paste(indicators_tree[-nber_indicators], collapse = "+")))
  
  selected_fit_glm <- glm(xformula, data = xdata, family = family)
  
  selected_p_value <- vect_p_values[index_best_tree]
  
  ######### outer layer of the permutation 
  
  if(family == "binomial")
  {
    Vect_PY1 <- fit_null$fitted
    Yb_list <- lapply(1:BB, function(j) 
    {
      return(sapply(Vect_PY1, function(xp) {
        return(sample(c(1, 0), 1, replace = TRUE, prob = c(xp, 1-xp)))
      }))
    })
  }
  else
  {
    df <- length(hat_gamma0)
    hat_sd2e <- sqrt(sum(fit_null$residuals^2)/(n-df))
    
    Yb_list <- list()
    Yb_list <- lapply(1:BB, function(j)
    {
      return(fit_null$fitted.values + rnorm(n, mean = 0, sd = hat_sd2e))
    })
  }
  
  
  ##	Parallel job
  numWorkers <- args.parallel$numWorkers
  cat("\n ncores = ", numWorkers, " for  outer layer bootstrap !\n")
  
  List_xdatas <- lapply(Yb_list, function(Ybb) return(data.frame(Yb = Ybb, xdata)))
  
  MaxTreeList <- mclapply(List_xdatas, wrapper, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  
  
  size_vect <- sapply(MaxTreeList, function(xxtree) {return(sum(xxtree$frame$var == "<leaf>"))})
  if (min(size_vect) <= m)
  {
    supp_size <- which(size_vect <= m)
    MaxTreeList[supp_size] <- NULL
    size_vect <- size_vect[-supp_size]
    Yb_list[supp_size] <- NULL
  }
  
  BB_adj <- length(size_vect)  
  
  
  List_xTrees_xDatas <- list()
  List_xTrees_xDatas <-	lapply(1:BB_adj, function(j)
  {
    return(list(MaxTreeList[[j]], data.frame(Yb = Yb_list[[j]], xdata)))
  })
  
  
  cat("\n ncores = ", numWorkers, " for nested trees in the outer layer !\n")
  
  List_Diff_deviancesBB <- mclapply(List_xTrees_xDatas, wrapper2, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  Diff_deviancesBB <- matrix(unlist(List_Diff_deviancesBB), nrow = m - 1, byrow = FALSE)
  mat_P_value <- t(sapply(1:(m-1), function(j){sapply(Diff_deviancesBB[j,], function(u){
    mean(Diff_deviances[j,] > u)})}))
  samp_p_value <- apply(mat_P_value,2,min)
  adj_p_value <- mean(samp_p_value < selected_p_value)
  
  ######## end of outer level of the permutation
  
  Tree_Selected <- (adj_p_value < LEVEL)
  
  selected_fit_glm_tree <- list(fit_glm = selected_fit_glm, tree = selected_tree, p.value = selected_p_value, adj_p.value = adj_p_value, Tree_Selected = Tree_Selected)
  
  if(! Tree_Selected) 
  {
    selected_fit_glm <- glm(as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "))), data = xdata, family = family)
  }
  time2 <- Sys.time()
  Timediff <- difftime(time2, time1)
  
  return(list(selected_model = selected_fit_glm_tree, fit_glm = selected_fit_glm, Timediff = Timediff, comp_p_values = vect_p_values, Badj = B_adj, BBadj = BB_adj))
}

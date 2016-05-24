p.val.tree <- function(xtree, xdata, Y.name, X.names, G.names, B = 10, args.rpart = list(minbucket = 40, maxdepth = 10, cp = 0), epsi = 1e-3, iterMax = 5, iterMin = 3, family = "binomial", LB = FALSE, args.parallel = list(numWorkers = 1), index = 4, verbose = TRUE)
{ 
  if(!inherits(xtree, 'rpart')) stop('xtree have to be an rpart object!')
  if (index <= 1) stop ("The test procedure is not available for a root tree")
  time1 <- Sys.time()
  ##	Fit null model
  fit_null <- glm(as.formula(paste(Y.name, " ~ ", paste(X.names, collapse = "+"))), data = xdata, family = family)
  
  ##	Parameter of X covariates
  hat_gamma0 <- fit_null$coef
  
  n <- nrow(xdata)
  
  ##	For product
  product <- ifelse(length(X.names) == 1, "*", "%*%")
  
  ##	Workers job
  MaxTreeList = list()
  
  
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
  }
  else
  {
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
  cat("\n ncores = ", numWorkers, " for bootstrap !\n")
  
  List_xdatas <- lapply(Yb_list, function(Ybb) return(data.frame(Yb = Ybb, xdata)))
  
  MaxTreeList <- mclapply(List_xdatas, wrapper, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  
  
  size_vect <- sapply(MaxTreeList, function(xxtree) {return(sum(xxtree$frame$var == "<leaf>"))})
  if (min(size_vect) <= index)
  {
    supp_size <- which(size_vect <= index)
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
  
  
  cat("\n ncores = ", numWorkers, " for nested trees !\n")
  
  List_Diff_deviances <- mclapply(List_xTrees_xDatas, wrapper2, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  
  ##	Observed statistics
  obs_nested_trees <- nested.trees(xtree, xdata, Y.name, X.names, MaxTreeSize = m, family = family, verbose = verbose)
  obs_diff_deviances <- obs_nested_trees$diff_deviances
  Diff_deviances <- matrix(unlist(List_Diff_deviances), nrow = m - 1, byrow = FALSE)
  
  ##  Compute p.values by Wang et al. procedure using the moment estimator
  mu <- apply(Diff_deviances,1,mean)
  sigma2 <- apply(Diff_deviances,1,var) 
  r <- 4*mu/sigma2
  bn <- r*mu/2
  trans_dev <- r*Diff_deviances/2
  trans_obs_dev <- r*obs_diff_deviances/2
  p_value_selected <- 1 - pchisq(trans_obs_dev[index-1], df = bn[index-1])
  
  time2 <- Sys.time()
  Timediff <- difftime(time2, time1)
  
  return(list(p.value = p_value_selected, Timediff = Timediff,Badj = B_adj))
  
}

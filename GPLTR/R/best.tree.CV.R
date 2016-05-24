## Cross validation
best.tree.CV <- function(xtree,xdata,Y.name, X.names, G.names, family = 'binomial', args.rpart = list(cp=0,minbucket=20,maxdepth=10), epsi = 1e-3, iterMax = 5, iterMin = 3, ncv = 10, verbose = TRUE)
{ 
  time1 <- Sys.time()
  if(!inherits(xtree, 'rpart')) stop('xtree must be an rpart object')
  tree_size1 <- sum(xtree$frame$var == "<leaf>")
  Ind_data <- 1:nrow(xdata)
  CV_ERRORS <- list()
  CV_ERROR0 <- list()
  for(l in 1:ncv)
  {
    if(l!=ncv)
    {
      Ind_test <- sample(Ind_data,size=round(nrow(xdata)/ncv))
      Ind_data <- Ind_data[!(Ind_data%in%Ind_test)]
    }
    else
    {
      Ind_test <- Ind_data
    }
    data_learn <- xdata[-Ind_test,]
    data_test <- xdata[Ind_test,]
    fit_pltr <- pltr.glm(data_learn, Y.name, X.names, G.names, family = family, args.rpart = args.rpart, epsi = epsi, iterMax = iterMax, iterMin = iterMin, verbose = verbose)
    fit0 <- glm(as.formula(paste(Y.name, "~ ", paste(X.names, collapse=" + "))), data = data_learn, family = family)
    pred0 <- predict.glm(fit0, newdata = data_test, type = "response")
    predict_glm0 <- as.numeric(pred0>0.5)
    CV_ERROR0[[l]] <- mean(data_test[Y.name] != predict_glm0)
    
    tree_size <- sum(fit_pltr$tree$frame$var == "<leaf>")
    if (tree_size > 1)
    {    
      MaxTreeSize <- min(tree_size, 10, tree_size1)
      nested_trees <- nested.trees(fit_pltr$tree, data_learn, Y.name, X.names, MaxTreeSize = MaxTreeSize, family = family, verbose = verbose)
      
      fits_glm <-  lapply(nested_trees$leaves, function(uu)
      {
        xxtree <- snip.rpart(fit_pltr$tree, toss = uu)
        xfit <- tree2glm(xxtree, data_learn, Y.name, X.names, family = family)
        return(xfit)  
      })
      
      predict_glm <- lapply(fits_glm,function(u)
      {
        pred <- predict.glm(u, newdata = data_test, type = "response")
        return(as.numeric(pred>0.5))
      })
      
      CV_ERRORS[[l]] <- sapply(predict_glm,function(u)
      { 
        return(mean(data_test[Y.name]!=u))
      })
    }
    else
    {
      CV_ERRORS[[l]] <- NULL
      CV_ERROR0[[l]] <- NULL
    }
    
  }
  CV_ERROR0 <- CV_ERROR0[!sapply(CV_ERROR0,is.null)]
  Vect_CV_ERROR0 <- mean(unlist(CV_ERROR0))
  CV_ERRORS <- CV_ERRORS[!sapply(CV_ERRORS,is.null)]
  dim_CV <- min(sapply(CV_ERRORS, length))
  CV_ERRORS <- lapply(CV_ERRORS, function(u) u[1:dim_CV])
  vect_CV_ERRORS <- apply(matrix(unlist(CV_ERRORS),ncol= length(CV_ERRORS),byrow=F),1,mean)
  CV_ERRORS <- c(Vect_CV_ERROR0,vect_CV_ERRORS)
  CV_ERROR <- min(CV_ERRORS)
  best_index_CV <- which.min(CV_ERRORS)
  nested_tree <- nested.trees(xtree, xdata, Y.name, X.names, family = family, verbose = verbose)
  nested_trees_leaves <- append(nested_tree$leaves, c(1), after = 0)
  best_tree_CV <- snip.rpart(xtree, toss = nested_trees_leaves[[best_index_CV]])
  if(best_index_CV > 1)
  {
    indicators_tree <- sapply(tree2indicators(best_tree_CV), function(u) return(paste("as.integer(", u, ")")))
    
    nber_indicators <- length(indicators_tree)
    
    xformula <- as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "), "+", paste(indicators_tree[-nber_indicators], collapse = "+")))
    
    lm_CV <- glm(xformula, data = xdata, family = family)
  }else{
    lm_CV <- glm(as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "))), data = xdata, family = family)
  }
  time2 <- Sys.time()
  Timediff <- difftime(time2, time1)
  return(list(best_index = best_index_CV, tree = best_tree_CV, fit_glm = lm_CV , CV_ERRORS = list(CV_ERROR,CV_ERRORS), Timediff = Timediff))
}

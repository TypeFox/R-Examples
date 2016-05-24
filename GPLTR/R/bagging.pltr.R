bagging.pltr <- function(xdata, Y.name, X.names, G.names, family = "binomial", args.rpart = list(minbucket = 20), epsi = 0.001, iterMax = 5, iterMin = 3, LB = FALSE, args.parallel = list(numWorkers = 1), Bag = 20, Pred_Data = data.frame(), verbose = TRUE, doprune = FALSE, thresshold = seq(0, 1, by = 0.1))
{
  time1 <- Sys.time()
  n = nrow(xdata)
  Ind_data <- 1:n
  IND_SAMP <- lapply(1:Bag, function(u) sample.int(n, size = n, replace = TRUE))
  IND_OOB  <- lapply(IND_SAMP, function(v) Ind_data[!(Ind_data %in% v)])
  
  wrapper <- function(xdata_bag)
  {
    pltr_lm_b = pltr.glm(xdata_bag, Y.name = Y.name, X.names = X.names, G.names = G.names, family = family, args.rpart = args.rpart, epsi = epsi, iterMax = iterMax, iterMin = iterMin, verbose = verbose)
    
    return(pltr_lm_b$tree)
  }
  numWorkers = args.parallel$numWorkers
  cat("\n ncores = ", numWorkers, " for bagging trees !\n")
  List_xdatas <- lapply(IND_SAMP, function(w) return(xdata[w,]))
  
  MaxTreeList <- mclapply(List_xdatas, wrapper, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
  
  List_xTrees_xDatas <-  lapply(1:Bag, function(j)
  {
    return(list(MaxTreeList[[j]], List_xdatas[[j]]))
  })
  if(doprune ){
    cat("\n ncores = ", numWorkers, " for bagging selected trees !\n")
    wrapper2 <- function(list_xtree_xdata)
    {
      resultBICAIC = best.tree.BIC.AIC(xtree = list_xtree_xdata[[1]], xdata = list_xtree_xdata[[2]], Y.name, X.names, family, verbose = verbose)
      return(list(resultBICAIC$tree$BIC, resultBICAIC$fit_glm$BIC))
    }
    tree_model_BIC <- mclapply(List_xTrees_xDatas, wrapper2, mc.cores = getOption("mc.cores", numWorkers), mc.preschedule = LB, mc.silent = TRUE)
    LIST_tree_BIC_Bag <- lapply(tree_model_BIC,function(mod) return(mod[[1]]))
    LIST_glm_Bag <- lapply(tree_model_BIC,function(modd){
      moda <- modd[[2]]
      moda$data <- NULL
      moda$model <- NULL
      moda$y <- NULL
      moda$fitted.values <- NULL
      moda$linear.predictors <- NULL
      moda$prior.weights <- NULL
      moda$weights <- NULL
      moda$qr$qr <- NULL
      return(moda)})
      ImpVarlist = lapply(LIST_tree_BIC_Bag, function(ww){
      ww$variable.importance
    })
    Impvarmat <- unlist(ImpVarlist)
    TREES <- LIST_tree_BIC_Bag
  }else{
    LIST_glm_Bag <- lapply( List_xTrees_xDatas, function(u){
      moda <- tree2glm(xtree = u[[1]], xdata = u[[2]], Y.name, X.names, family)
      moda$data <- NULL
      moda$model <- NULL
      moda$y <- NULL
      moda$fitted.values <- NULL
      moda$linear.predictors <- NULL
      moda$prior.weights <- NULL
      moda$weights <- NULL
      return(moda)
    })

    ImpVarlist = lapply(MaxTreeList, function(ww){
      ww$variable.importance
    })
    Impvarmat <- unlist(ImpVarlist)
    TREES <- MaxTreeList
  }
  ImpVar <- tapply(Impvarmat, names(Impvarmat), sum)
  ImpVar <- ImpVar/sum(ImpVar)*100
  ImpVar <- round(sort(ImpVar, decreasing = TRUE), 2)
  
  List_xdatas_OOB <- lapply(IND_OOB, function(w) return(xdata[w,]))
  List_xdatas_Glm_OBB <- lapply(1:Bag, function(j)
  {
    return(list(List_xdatas_OOB[[j]], LIST_glm_Bag[[j]]))
  })
  predict_glm_OOB_PBP <- lapply(List_xdatas_Glm_OBB,function(uw)
  {
    pred = predict.glm(uw[[2]], newdata = uw[[1]], type = "response")
    return(sapply(thresshold, function(vv) as.numeric(pred >= vv)))
  })
  OOB_ERRORS_PBP <- sapply(1:Bag,function(uuu)
  { 
    return(apply(predict_glm_OOB_PBP[[uuu]], 2, function(wz) mean(List_xdatas_OOB[[uuu]][Y.name] != wz)))
  })
  if(length(thresshold) > 1){
  OOB_ERROR_PBP <- apply(OOB_ERRORS_PBP, 1, mean)
  } else OOB_ERROR_PBP <- mean(OOB_ERRORS_PBP)
  
    
  ## Compute the OOB error for each OOB individual
  UNIQUE_IND_OOB <- sort(unique(unlist(IND_OOB)))
  LOST <- matrix(rep(0,length(UNIQUE_IND_OOB)*length(thresshold)), ncol = length(thresshold))
  j <- 0
  for(i in UNIQUE_IND_OOB){
    j <- j+1
    poslist <- sapply(IND_OOB, function(uu) is.element(i,uu))
    IND_OOBi <- IND_OOB[poslist]
    PRED_OOBi <- predict_glm_OOB_PBP[poslist]
    posi <- sapply(IND_OOBi, function(w) which(w == i))

    vec_PREDi <- sapply(seq(length(posi)), function(ww) PRED_OOBi[[ww]][posi[ww],])
    if(length(thresshold) >1){
    PREDi <- apply(vec_PREDi, 1, function(zz) as.numeric(mean(zz) > 0.5))
    } else PREDi <- as.numeric(mean(vec_PREDi) > 0.5)
    LOST[j, ] <- sapply(PREDi, function(ss) ss!= t(xdata[Y.name])[i]) 
  }
  EOOB <- apply(LOST, 2, mean)
  names(EOOB) <- paste('CUT', 1:length(thresshold), sep = '')
  
  if(nrow(Pred_Data) != 0)
  {
    predict_glm2 <- lapply(LIST_glm_Bag, function(uw)
    {
      pred = predict.glm(uw, newdata = Pred_Data, type = "response")
      return(sapply(thresshold, function(wz) as.numeric(pred > wz)))
    })
    
    PRED_IND <- list()
    for(jj in seq(length(thresshold))){
      PRED_IND[[jj]] <- sapply(1:Bag, function(ww) predict_glm2[[ww]][,jj])
    }
    
    FINAL_PRED_IND <- lapply(PRED_IND, function(www) apply(www, 1, function(zzz) as.numeric(mean(zzz) > 0.5)))
    PRED_ERROR <- sapply(FINAL_PRED_IND, function(uuu) mean( uuu != Pred_Data[Y.name]))    
  }
  time2 <- Sys.time()
  Timediff <- difftime(time2, time1)
  return(list(IND_OOB = IND_OOB, EOOB = EOOB, OOB_ERRORS_PBP = OOB_ERRORS_PBP, OOB_ERROR_PBP = OOB_ERROR_PBP,
              Tree_BAG = TREES, Glm_BAG = LIST_glm_Bag, LOST = LOST, TEST = (if(length(Pred_Data) != 0)
              list(PRED_ERROR = PRED_ERROR, PRED_IND = PRED_IND, FINAL_PRED_IND = FINAL_PRED_IND) else NULL), 
              Var_IMP = ImpVar, Timediff = Timediff, CUT = thresshold))  
}
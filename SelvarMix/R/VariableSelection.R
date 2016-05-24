##################################################################
##          VariableSelection.R
##################################################################
VariableSelection<-
  function(data,
           nbCluster,
           models,
           criterion,
           OrderVariable,
           hybrid.size,
           supervised,
           knownlabels,
           nbCores)
  {
    data <- as.matrix(data)
    nbCluster <- as.integer(nbCluster)
    criterion <- as.character(criterion)
    
    
    listModels.size <- length(models["listModels"])
    nbCluster.size <- length(nbCluster)
    OutputVector.size <- listModels.size*nbCluster.size
    
    wrapper.selectVar <- function(arg)
    {
      myModel <-  mixmodGaussianModel(listModels=models["listModels"][arg[2]])
      orderVar <- rep(NA, ncol(data))
      
      if(length(nbCluster)==1)
        orderVar <- OrderVariable 
      else
        orderVar <- OrderVariable[arg[1],]
      
      my.list <- rcppSelectS(data, 
                             orderVar, 
                             nbCluster[arg[1]], 
                             myModel, 
                             hybrid.size, 
                             criterion,
                             knownlabels,
                             supervised)
      
      ResSelectVar <- list()
      ResSelectVar$S <- my.list$S
      ResSelectVar$nbCluster <- my.list$nbCluster
      ResSelectVar$criterion <- my.list$criterion
      ResSelectVar$model <- my.list$model 
      ResSelectVar$criterionValue  <- my.list$criterionValue
      ResSelectVar$partition <- my.list$partition
      ResSelectVar$proba <- my.list$proba
      OrderAux <- setdiff(OrderVariable, ResSelectVar$S)
      ResSelectVar$W <- rcppSelectW(data, OrderAux, ResSelectVar$S, hybrid.size)
      return(ResSelectVar)
    }
   
    if(OutputVector.size < nbCores)
      nbCores <- OutputVector.size
    
    arg.grid <- matrix(0, OutputVector.size, 2)  
    arg.grid <- as.matrix(expand.grid(1:nbCluster.size, 1:listModels.size))
    ## si on est sous windows
    if(Sys.info()["sysname"] == "Windows")
    {
      cl <- makeCluster(nbCores)
      common.objects <- c("data", 
                          "OrderVariable", 
                          "nbCluster",
                          "models",
                          "hybrid.size", 
                          "criterion",
                          "supervised",
                          "knownlabels",
                          "mixmodLearn",
                          "mixmodCluster",
                          "mixmodStrategy")
      #clusterEvalQ(cl, require(Rmixmod))
      clusterExport(cl=cl, varlist = common.objects, envir = environment())
      junk <- parApply(cl = cl,  
                        X = arg.grid,
                        MARGIN = 1,
                        FUN = wrapper.selectVar)
      stopCluster(cl)
      
    }
    else
      junk <- mclapply(X = as.list(data.frame(t(arg.grid))),
                       FUN = wrapper.selectVar,
                       mc.cores = nbCores,
                       mc.silent = FALSE,
                       mc.preschedule = TRUE,
                       mc.cleanup = TRUE)
    
    nb.fails <- 0 
    for(idx in 1:OutputVector.size)
      if(class(junk[[idx]]) == "try-error")
        nb.fails <- nb.fails + 1
    
    VariableSelectRes <-  vector(length = (OutputVector.size - nb.fails), mode ="list")
    
    idx <- 1
    for(ll in 1:OutputVector.size)
      if(class(junk[[ll]])!="try-error")
      {
        VariableSelectRes[[idx]]$S <- sort(junk[[ll]][["S"]])
        VariableSelectRes[[idx]]$W <- sort(junk[[ll]][["W"]])
        VariableSelectRes[[idx]]$U <- setdiff(1:dim(data)[2], union(junk[[ll]][["S"]], junk[[ll]][["W"]]))
        VariableSelectRes[[idx]]$criterionValue <- junk[[ll]][["criterionValue"]]
        VariableSelectRes[[idx]]$criterion <- junk[[ll]][["criterion"]] 
        VariableSelectRes[[idx]]$model <- junk[[ll]][["model"]]
        VariableSelectRes[[idx]]$nbCluster <- junk[[ll]][["nbCluster"]]
        VariableSelectRes[[idx]]$partition <- junk[[ll]][["partition"]]
        VariableSelectRes[[idx]]$proba <- junk[[ll]][["proba"]]
        idx <- idx + 1
      }
    
    return(VariableSelectRes)  
  }



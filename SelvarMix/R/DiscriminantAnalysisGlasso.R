DiscriminantAnalysisGlasso <- function(data, 
                                       nbCluster, 
                                       lambda, 
                                       rho,
                                       knownlabels,
                                       nbCores)
  {
  data <- as.matrix(data)
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  nbCluster <- as.integer(nbCluster)
  
  if((length(lambda)*length(rho)) < nbCores)
    nbCores <- (length(lambda)*length(rho))
  #print(c(" ... nbCores = ... ", nbCores))
  wrapper.DiscriminantAnalysisGlasso <- function(prm)
  {
    result <- rcppDiscriminantAnalysisGlasso(data, knownlabels, nbCluster, prm[1], prm[2])
    return(result)
  }
  
  pen.grid <- matrix(0, (length(lambda)*length(rho)), 2)  
  pen.grid <- as.matrix(expand.grid(lambda, rho))
#   pen.grid.list <- list(); colnames(pen.grid) <- NULL
#   pen.grid.list <- as.list(data.frame(t(pen.grid)))
#   
  ## si on est sous windows
  if(Sys.info()["sysname"] == "Windows")
  {
    cl <- makeCluster(nbCores)
    common.objects <- c("data", "nbCluster", "knownlabels","glasso") 
    #clusterEvalQ(cl, require(glasso))
    clusterExport(cl=cl, varlist = common.objects, envir = environment())
    parallel.varrole <-  parApply(cl = cl, 
                                   X = pen.grid,
                                   MARGIN =  1,
                                   FUN = wrapper.DiscriminantAnalysisGlasso)
    stopCluster(cl)
  }
  else
  parallel.varrole <-  mclapply(X = as.list(data.frame(t(pen.grid))),
                                FUN = wrapper.DiscriminantAnalysisGlasso,
                                mc.cores = nbCores,
                                mc.preschedule = TRUE,
                                mc.cleanup = TRUE )
 
  var.role <- matrix(0,(length(lambda)*length(rho)), p)
  for(j in 1:nrow(var.role))
    if(class(parallel.varrole[[j]])!="try-error")
      var.role[j,] <- parallel.varrole[[j]]   
  
  return(var.role)
}
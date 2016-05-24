Analyze_oneAE <- function(ae,
                          drug,
                          maxit,
                          alpha,
                          nbinit)
{
  if (sum(ae)>1){
    datau <- uniquecombs(cbind(ae,drug[,-which(colSums(drug[which(ae==1),])==0)]))
    w <- table(attr(datau,"index"))
    xx <- datau[,-1]
    y <- datau[,1]
    
    drugs.names <- colnames(drug[,-which(colSums(drug[which(ae==1),])==0)])
    
    if (ncol(datau) > 12){
      cat("Metropolis-Hastings is used since many drugs are considered\n")
      
      maxitlist <- as.list(rep(maxit,nbinit))
      nb.cpus <- min(detectCores(all.tests = FALSE, logical = FALSE) , nbinit)
      if(Sys.info()["sysname"] == "Windows")
      {
        cl <- makeCluster(nb.cpus)
        common.objects <- c("y","xx","w", "drugs.names", "alpha") 
        clusterExport(cl=cl, varlist = common.objects, envir = environment())
        results <- parLapply(cl = cl,
                             X = maxitlist,
                             fun = MHLogisticw,
                             y = y,
                             xx = xx,
                             w = w,
                             drugs.names = drugs.names,
                             alpha = alpha)
        #stop("Parallelisation is not available for windows")
      }
      else
      {
        results <- mclapply(X = maxitlist,
                            FUN = MHLogisticw,
                            y = y,
                            xx = xx,
                            w = w,
                            drugs.names = drugs.names,
                            alpha = alpha,
                            mc.cores = nb.cpus,
                            mc.preschedule = TRUE,
                            mc.cleanup = TRUE
        )
       
      }  
      idx <- NA
      bic <- rep(NA, length(results))
      for (k in 1:length(results)) bic[k] <- results[[k]]$best.model.bic
      results <- list(best=results[[which.max(bic)]], all=results)
      results$signals <- FindSignals(results)
    }else{
      cat("Exhaustive approach is used since few drugs are considered\n")
      results <- list(best=ExhaustiveLogisticw(y = y, x = as.matrix(xx), w=w, drugs.names = drugs.names))
      results$signals <- FindSignals(results)
    }
  }else{
    cat("The adverse event is notified only one time, so the estimation is not doable!\n")
    results <- list(best=list(model=NULL))
  }
  
  return(results)  
}



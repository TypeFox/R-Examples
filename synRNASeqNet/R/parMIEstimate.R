parMIEstimate <-
function(counts, method = c("ML", "MM", "Bayes", "CS",
                                             "Shrink", "KD", "KNN"),
                          unit = c("bit", "ban", "nat"), nchips, priorHyperParam =
                            c("Jeffreys", "BLUnif", "Perks", "MiniMax"),
                          shrinkageTarget, k = 3, tfList = NULL, boot = F){  
  method <- match.arg(method)
  unit <- match.arg(unit)
  
  if(method == "ML" || method == "MM" || method == "Bayes" || method == "CS" ||
       method == "Shrink"){
    g <- nrow(counts)
        
    if(missing(nchips)) nchips <- detectCores() - 1
    cl <- makeCluster(nchips)
    
    dataToExport <- "counts"
    
    
    if(missing(shrinkageTarget)) shrinkageTarget <- NULL
    
    if(method == "ML") funcToExport <- c("thetaML", "entropyML")
    if(method == "MM") funcToExport <- c("thetaML", "entropyML", "entropyMM")
    if(method == "Bayes"){
      funcToExport <- c("beta_k", "thetaBayes", "entropyBayes")
      priorHyperParam <- match.arg(priorHyperParam)
    }
    if(method == "CS") funcToExport <- c("thetaML", "thetaGT", "entropyCS")
    if(method == "Shrink") funcToExport <- c("thetaML", "shrinkageIntensity",
                                             "thetaShrink", "entropyShrink")
    
    toExport <- c(dataToExport, funcToExport)
    clusterExport(cl, toExport)
    
    if(is.null(tfList)){
      tmpEntropy <- parSapply(cl, 1:((g^2 + g)/2), parEntropyEstimate, method = method,
                              unit = unit, priorHyperParam = priorHyperParam,
                              shrinkageTarget, boot = boot)
            
      hatEntropy <- diag(0, g)
      hatEntropy[upper.tri(hatEntropy, diag = T)] <- tmpEntropy
      hatEntropy[lower.tri(hatEntropy)] <- t(hatEntropy)[lower.tri(t(hatEntropy))]
      
      MIEst <- function(i, ii)
        hatEntropy[i, i] + hatEntropy[ii, ii] - hatEntropy[i, ii]
      
      ans <- diag(0, g)
      for(i in 1:g)
        for(ii in i:g)
          ans[i, ii] <- MIEst(i, ii)
      ans[lower.tri(ans)] <- t(ans)[lower.tri(t(ans))]
      rownames(ans) <- colnames(ans) <- rownames(counts)
    }else{
      sstart <- startId(nrow(counts))
      eend <- diagId(nrow(counts))
      
      tfIdx <- which(rownames(counts) %in% tfList)
      tasks <- NULL
      for(tf in tfIdx){
        tasks <- c(tasks, sstart[tf]:eend[tf])
        if(tf != g) for(i in 0:(g-tf-1)) tasks <- c(tasks, tasks[length(tasks)] + tf + i)
      }
      
      tmpEntropy <- parSapply(cl, c(tasks, eend), parEntropyEstimate, method = method,
                              unit = unit, priorHyperParam = priorHyperParam,
                              shrinkageTarget, boot = boot)
      
      hatEntropy <- matrix(tmpEntropy[1:(length(tmpEntropy)-g)], nrow = length(tfList), ncol = g, byrow = T)
      tmpEntropy <- tmpEntropy[(length(tmpEntropy)-g+1):length(tmpEntropy)]
      
      MIEst <- function(i, ii)
        tmpEntropy[tfIdx[i]] + tmpEntropy[ii] - hatEntropy[i, ii]
      
      ans <- matrix(0, nrow = length(tfList), ncol = g)
      for(i in 1:length(tfIdx))
        for(ii in 1:g)
          ans[i, ii] <- MIEst(i, ii)
      rownames(ans) <- tfList
      colnames(ans) <- rownames(counts)
    }
    
    stopCluster(cl)
  } else
    if(method == "KD"){
      g <- nrow(counts)
      
      if(missing(nchips)) nchips <- detectCores() - 1
      cl <- makeCluster(nchips)
      
      dataToExport <- "counts"
      
      funcToExport <- c("parMIKD")
      toExport <- c(dataToExport, funcToExport)
      clusterExport(cl, toExport)
      
      if(is.null(tfList)){
        tmpMI <- parSapply(cl, 1:((g^2 + g)/2), parMIKD)#, unit = unit)
        
        ans <- diag(0, g)
        ans[upper.tri(ans, diag = T)] <- tmpMI
        ans[lower.tri(ans)] <- t(ans)[lower.tri(t(ans))]
        rownames(ans) <- colnames(ans) <- rownames(counts)
      }else{
        sstart <- startId(nrow(counts))
        eend <- diagId(nrow(counts))
        
        tfIdx <- which(rownames(counts) %in% tfList)
        tasks <- NULL
        for(tf in tfIdx){
          tasks <- c(tasks, sstart[tf]:eend[tf])
          if(tf != g) for(i in 0:(g-tf-1)) tasks <- c(tasks, tasks[length(tasks)] + tf + i)
        }
        
        tmpMI <- parSapply(cl, tasks, parMIKD)#, unit = unit)
        
        ans <- matrix(tmpMI, nrow = length(tfList), ncol = g, byrow = T)
        rownames(ans) <- tfList
        colnames(ans) <- rownames(counts)
      }
      
      stopCluster(cl)
    } else
      if(method == "KNN"){
        if(is.null(tfList)) ans <- knnmi.cross(counts, counts, k = k, noise = 0) else
        #if(is.null(tfList)) ans <- knnmi.all(counts, k = k, noise = 0) else
          ans <- knnmi.cross(counts[tfList, ], counts, k = k, noise = 0)
      } else
        stop("Unknown Mutual Information Estimate Method")
  
  return(ans)
}

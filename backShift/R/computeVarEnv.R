computeVarEnv <- function(B, DeltaList, base.setting.index = 1, verbose = FALSE){
  
  if(is.null(base.setting.index)) base.setting.index <- 1
  
  if(verbose){
    cat("Estimating intervention variance... \n")
  }
  
  G <- length(DeltaList)
  DiagList <- vector("list", G)
  for(i in 1:G){
    DiagList[[i]] <- B %*% DeltaList[[i]] %*% t(B)
  }
  
  p <- nrow(B)
  
  relIntVars <- matrix(0, G, p)
  base.setting <- DiagList[[base.setting.index]]
  for(i in 1:G){
    if(i != base.setting.index){
      relIntVars[i,] <- (G-1)/G * diag((DiagList[[i]] - base.setting))
    }
  }
  if(verbose){
    cat("Variance of interventions are estimated relative to setting", base.setting.index, "\n")
    cat("Estimating intervention variance...done! \n")
  }
  
  relIntVars
}
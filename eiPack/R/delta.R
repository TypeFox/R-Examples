delta <- function(object, cols){
  se <- matrix(NA, nrow(object[["se"]]), length(cols))
  rownames(se) <- rownames(object[["se"]])
  colnames(se) <- cols
  cov.array <- array(dim=c(nrow(se), nrow(se), ncol(se)))
  rownames(cov.array) <- colnames(cov.array) <- rownames(se)
  dimnames(cov.array)[[3]] <- cols
  for(i in 1:length(cols)){
    cov.array[,,i] <- object[["cov.matrices"]][cols[i]][[1]]
  }
                           
  delta.nums <- paste("x",1:length(cols), sep="")
  delta.denom <-  paste(delta.nums[1],"+", delta.nums[2], sep="")
  if(length(delta.nums)>2){
    for(i in 3:length(delta.nums)){
      delta.denom <- paste(delta.denom, "+", delta.nums[i], sep="")
    }
  }
  
  delta.formula <- list()
  for(i in 1:length(delta.nums)){
    delta.formula[[i]] <- as.formula(paste("~ ", delta.nums[i], "/(",
                                           delta.denom,")", sep=""))
  }
  

  for(i in 1:nrow(se)){
      se[i,] <- msm::deltamethod(delta.formula,
                                  m=object[["coefficients"]][i,cols],
                                  cov=diag(cov.array[i,i,cols]))
    }
  se
}

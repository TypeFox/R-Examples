residuals.ivmodel <- function(object, ...){
  ivmodel <- object
  temp <- matrix(NA, nrow=2, ncol=0)
  if(!is.null(ivmodel$kClass)){
    temp<-rbind(as.numeric(rownames(ivmodel$kClass$ci)), c(ivmodel$kClass$point.est))
	colnames(temp) <- paste("k=", temp[1,], sep="")
    colnames(temp)[temp[1,]==0] <- "OLS"
    colnames(temp)[temp[1,]==1] <- "TSLS"
  }
  if(!is.null(ivmodel$LIML)){  
    LIML <- c(ivmodel$LIML$k, ivmodel$LIML$point.est)
	temp <- cbind(temp, LIML)
  }
  if(!is.null(ivmodel$Fuller)){  
    Fuller <- c(ivmodel$Fuller$k, ivmodel$Fuller$point.est)
	temp <- cbind(temp, Fuller)
  }

  result <- matrix(NA, ncol=ncol(temp), nrow=ivmodel$n)
  colnames(result) <- colnames(temp)
  for(i in 1:ncol(temp)){
    result[, i] <- ivmodel$Y - ivmodel$Yadj - temp[2, i]*(ivmodel$D-ivmodel$Dadj)
  }
  return(result)
}
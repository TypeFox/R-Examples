coef.ivmodel<-function(object, ...){
  ivmodel<-object
  coefmat <- matrix(NA, ncol=5, nrow=0)
  colnames(coefmat) <- c("k", "Estimate", "Std. Error", "t value", "Pr(>|t|)")
  if(!is.null(ivmodel$kClass)){
    temp<-cbind(as.numeric(rownames(ivmodel$kClass$ci)), ivmodel$kClass$point.est, 
            	ivmodel$kClass$std.err, ivmodel$kClass$test.stat, ivmodel$kClass$p.value)
	rownames(temp) <- rep("k-class", nrow(temp))
    rownames(temp)[temp[,1]==0] <- "OLS"
    rownames(temp)[temp[,1]==1] <- "TSLS"
    coefmat <- rbind(coefmat, temp)
  }
  if(!is.null(ivmodel$LIML)){
    temp<-cbind(ivmodel$LIML$k, ivmodel$LIML$point.est, ivmodel$LIML$std.err,
	            ivmodel$LIML$test.stat, ivmodel$LIML$p.value)
	rownames(temp) <- "LIML"
    coefmat <- rbind(coefmat, temp)
  }
  if(!is.null(ivmodel$Fuller)){
    temp<-cbind(ivmodel$Fuller$k, ivmodel$Fuller$point.est, ivmodel$Fuller$std.err, ivmodel$Fuller$test.stat, ivmodel$Fuller$p.value)
	rownames(temp) <- "Fuller"
    coefmat <- rbind(coefmat, temp)
  }
  
  coefmat<-coefmat[sort(coefmat[,1], index.return=T)$ix,]  
  return(coefmat)
}
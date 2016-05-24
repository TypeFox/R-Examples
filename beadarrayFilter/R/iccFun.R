iccFun <-
function(eSet,seSet,nSet,ProbeID,iccQuant,diffIcc=TRUE,keepData=TRUE)
{
    group <- c(1:dim(eSet)[2])
    fit1 <- MLM.beadarray(eSet, seSet, nSet, list(group), var.equal = TRUE,max.iteration = 20, method = "ML")
    est_var <-  fit1[, grep("sigma2", names(fit1))]
    tau <- fit1$tau2
print("Computing ICC for each array")
    iccAll <- sapply(1:dim(est_var)[2], function(i) tau/(tau + est_var[, i]))
###if iccQuant is a vector
print("Summarising the ICC at supplied quantile(s)")
if (diffIcc) 
{
###when diffIcc=TRUE, iccQuant must be a vector, give a warning here
if(length(iccQuant)<2) {stop("iccQuant must be a vector,  since diffIcc=TRUE. Specify diffIcc=FALSE")}
else if(length(iccQuant)>1)  {
    icc <- t(sapply(1:nrow(iccAll), function(x) quantile(iccAll[x,], iccQuant)))
   icc1 <- data.frame(ProbeID, icc)
   colnames(icc1) <- c( "ProbeID", paste("q", iccQuant,  sep = ""))
}
}
else{ 
###if iccQuant is a single value
icc <- sapply(1:nrow(iccAll), function(x) quantile(iccAll[x,], iccQuant))
icc1 <- data.frame(ProbeID, icc)
    colnames(icc1) <- c( "ProbeID", paste("q", iccQuant,  sep = ""))
}
 est_var1 <- data.frame( ProbeID, fit1[, grep("sigma2", names(fit1))])
    tau1 <- data.frame( ProbeID, fit1$tau2)
if (keepData)
       return(list(icc=icc1,withinvar=est_var1,betweenvar=tau1,iccAll =iccAll))
else  return(list(icc=icc1, withinvar=NULL,betweenvar=NULL,iccAll =NULL))
}

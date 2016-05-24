GeneratePvalues <-
function(data, B, addit=FALSE, covariable=NULL, family=binomial) {
  pvalors <- numeric(ncol(data)-1)
  pvalors <- sort(runPvalues(data, addit, covariable,family))
  pvalors <- rbind(pvalors, t(sapply(1:B, function(i) sort(runPermut(data, addit, covariable, family)))))
  rownames(pvalors) <- NULL  
  return(pvalors)
}

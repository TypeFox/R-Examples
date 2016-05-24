prop.multcomp <-
function(x,p,p.method="fdr") {
  obs <- integer(nrow(x))
  pval <- integer(nrow(x))
  for (i in 1:nrow(x)) {
    obs[i] <- x[i,1]/sum(x[i,])
    test <- binom.test(x[i,1],sum(x[i,]),p[i])
    pval[i] <- test$p.value
  }
  p.adj <- p.adjust(pval,method=p.method)
  comp <- data.frame("observed"=obs,"expected"=p,"p-value"=p.adj," "=.psignif(p.adj),stringsAsFactors=FALSE,check.names=FALSE)
  if (!is.null(rownames(x))) {rownames(comp) <- rownames(x)}
  dname <- paste(quote(x)," and ",quote(p),sep="")
  result <- list(method="exact binomial tests",data.name=dname,observed=obs,expected=p,p.adjust.method=p.method,p.value2=p.adj,p.value=comp)
  class(result) <- "RV.multcomp"
  return(result)
}

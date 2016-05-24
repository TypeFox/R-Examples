spearman.cor.multcomp <-
function(var1,var2,fact,alpha=0.05,nrep=1000){
  if (length(var1)!=length(var2)) {stop("'var1' and 'var2' lengths differ")}
  if (length(var1)!=length(fact)) {stop("'var1' and 'fact' lengths differ")}
  if (length(var2)!=length(fact)) {stop("'var2' and 'fact' lengths differ")}
  if (!is.factor(fact)) {fact <- as.factor(fact)}
  var1.2 <- var1[complete.cases(var1,var2,fact)]
  var2.2 <- var2[complete.cases(var1,var2,fact)]
  fact.2 <- droplevels(fact[complete.cases(var1,var2,fact)])
  dname <- paste(deparse(substitute(var1))," and ",deparse(substitute(var2))," by ",deparse(substitute(fact)),sep="")
  nlev <- nlevels(fact.2)
  tab <- data.frame(inf=integer(nlev),r=integer(nlev),sup=integer(nlev),row.names=levels(fact.2))
  cl <- 1-(alpha/nlev)
  for (i in 1:nlev) {
    ci <- spearman.ci(var1.2[as.numeric(fact.2)==i],var2.2[as.numeric(fact.2)==i],conf.level=cl,nrep=nrep)
    tab[i,] <- c(ci$conf.int["Inf"],ci$estimate,ci$conf.int["Sup"])
  }
  met <- paste("Comparison of ",nlev," Spearman's correlation coefficients",sep="")
  result <- list(method=met,data.name=dname,tab=tab,alpha=alpha,nrep=nrep)
  class(result) <- c("spearman.cor.multcomp","list")
  return(result)
}

print.spearman.cor.multcomp <- function(x,...) {
  cat("\n")
  cat(strwrap(x$method,prefix="\t"),sep="\n")
  cat("\n")
  cat("data: ",x$data.name,"\n")
  cat("Bonferroni-adjusted ",round(100*(1-x$alpha),1),"% confidence intervals\n")
  cat(x$nrep," replicates\n\n")
  print(x$tab,digits=4)
}

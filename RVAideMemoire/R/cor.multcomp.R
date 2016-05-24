cor.multcomp <-
function(var1,var2,fact,alpha=0.05,conf.level=0.95,theo=0,p.method="fdr"){
  if (length(var1)!=length(var2)) {stop("'var1' and 'var2' lengths differ")}
  if (length(var1)!=length(fact)) {stop("'var1' and 'fact' lengths differ")}
  if (length(var2)!=length(fact)) {stop("'var2' and 'fact' lengths differ")}
  if (theo<(-1) | theo>1) {stop("'theo' must be between -1 and 1")}
  if (!is.factor(fact)) {fact <- as.factor(fact)}
  nul <- as.numeric(row.names(table(c(which(is.na(var1)),which(is.na(var2))))))
  var1.2 <- if(length(nul)>0) {var1[-nul]} else {var1}
  var2.2 <- if(length(nul)>0) {var2[-nul]} else {var2}
  fact.2 <- if(length(nul)>0) {fact[-nul]} else {fact}
  dname <- paste(deparse(substitute(var1))," and ",deparse(substitute(var2))," by ",deparse(substitute(fact)),sep="")
  n <- integer(nlevels(fact.2))
  r <- integer(nlevels(fact.2))
  names(r) <- paste("coeff in group ",levels(fact),sep="")
  z <- integer(nlevels(fact.2))
  for (i in 1:nlevels(fact.2)) {
    n[i] <- length(var1.2[fact.2==levels(fact.2)[i]])
    r[i] <- as.numeric(cor.test(var1.2[fact.2==levels(fact.2)[i]],var2.2[fact.2==levels(fact.2)[i]])$estimate)
    z[i] <- 0.5*log((1+r[i])/(1-r[i]))
  }
  z.moy <- sum((n-3)*z)/sum(n-3)
  chi2.obs <- sum((n-3)*(z-z.moy)^2)
  names(chi2.obs) <- "X-squared"
  p <- pchisq(chi2.obs,nlevels(fact)-1,lower.tail=FALSE)
  met <- paste("Comparison of ",nlevels(fact)," Pearson's linear correlation coefficients",sep="")
  nval <- 0
  names(nval) <- "difference in coefficients"
  result <- list(method.test=met,data.name=dname,statistic=chi2.obs,parameter=c("df"=nlevels(fact)-1),p.value=p,
    alternative="two.sided",null.value=nval,estimate=r,alpha=alpha)
   if (p>alpha) {
    r.com <- (exp(2*z.moy)-1)/(exp(2*z.moy)+1)
    z.moy.inf <- z.moy-qnorm((1+conf.level)/2,0,1)/sqrt(sum(n)-3*nlevels(fact))
    z.moy.sup <- z.moy+qnorm((1+conf.level)/2,0,1)/sqrt(sum(n)-3*nlevels(fact))
    r.com.inf <- (exp(2*z.moy.inf)-1)/(exp(2*z.moy.inf)+1)
    r.com.sup <- (exp(2*z.moy.sup)-1)/(exp(2*z.moy.sup)+1)
    zeta <- 0.5*log((1+theo)/(1-theo))
    u.obs.com <- abs(z.moy-zeta)*sqrt(sum(n)-3*nlevels(fact))
    names(u.obs.com) <- "u"
    p.com <- min(pnorm(u.obs.com,0,1),pnorm(u.obs.com,0,1,lower.tail=FALSE))*2
    tab.com<-data.frame("inf"=r.com.inf,"r"=r.com,"sup"=r.com.sup,"theoretical"=theo,"U"=u.obs.com,"Pr(>|U|)"=p.com,
	" "=.psignif(p.com),stringsAsFactors=FALSE,check.names=FALSE)
    result$conf.level <- conf.level
    result$common.name <- paste("        Common correlation coefficient, ",100*conf.level,"% confidence interval\n",
	"          and equality to given value ",theo,sep="")
    result$common <- tab.com
  }
  if (p<alpha & nlevels(fact.2)>2) {
    fun.p <- function(i,j) {
	u <- abs(z[i]-z[j])/sqrt(1/(n[i]-3)+1/(n[j]-3))
	min(pnorm(u,0,1),pnorm(u,0,1,lower.tail=FALSE))*2
    }
    result$p.adjust.method <- p.method
    result$p.value.multcomp <- pairwise.table(fun.p,levels(fact.2),p.adjust.method=p.method)
  }
  class(result) <- "RVtest"
  return(result)
}

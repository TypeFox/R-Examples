cor.2comp <-
function(var1,var2,var3,var4,alpha=0.05,conf.level=0.95,theo=0){
  if (length(var1)!=length(var2)) {stop("'var1' and 'var2' lengths differ")}
  if (length(var3)!=length(var4)) {stop("'var3' and 'var4' lengths differ")}
  if (theo<(-1) | theo>1) {stop("'theo' must be between -1 and 1")}
  nul1 <- as.numeric(row.names(table(c(which(is.na(var1)),which(is.na(var2))))))
  var1.2 <- if(length(nul1)>0) {var1[-nul1]} else {var1}
  var2.2 <- if(length(nul1)>0) {var2[-nul1]} else {var2}
  nul2 <- as.numeric(row.names(table(c(which(is.na(var3)),which(is.na(var4))))))
  var3.2 <- if(length(nul2)>0) {var3[-nul2]} else {var3}
  var4.2 <- if(length(nul2)>0) {var4[-nul2]} else {var4}
  r1 <- as.numeric(cor.test(var1.2,var2.2,method="pearson")$estimate)
  r2 <- as.numeric(cor.test(var3.2,var4.2,method="pearson")$estimate)
  z1 <- 0.5*log((1+r1)/(1-r1))
  z2 <- 0.5*log((1+r2)/(1-r2))
  u.obs <- abs(z1-z2)/sqrt(1/(length(var1.2)-3)+1/(length(var3.2)-3))
  names(u.obs) <- "U"
  p <- 2*min(pnorm(u.obs,0,1),pnorm(u.obs,0,1,lower.tail=FALSE))
  met <- "Comparison of 2 Pearson's linear correlation coefficients"
  nval <- 0
  names(nval) <- "difference in coefficients"
  estimate <- c(r1,r2)
  names(estimate) <- paste("coeff in group ",1:2,sep="")
  result <- list(method.test=met,data.name="4 variables",statistic=u.obs,p.value=p,alternative="two.sided",
    null.value=nval,estimate=estimate,alpha=alpha)
  if (p>alpha){
    z.moy <- sum((length(var1.2)-3)*z1,(length(var3.2)-3)*z2)/sum(length(var1.2)-3,length(var3.2)-3)
    r.com <- (exp(2*z.moy)-1)/(exp(2*z.moy)+1)
    z.moy.inf <- z.moy-qnorm((1+conf.level)/2,0,1)/sqrt(sum(length(var1.2),length(var3.2))-6)
    z.moy.sup <- z.moy+qnorm((1+conf.level)/2,0,1)/sqrt(sum(length(var1.2),length(var3.2))-6)
    r.com.inf <- (exp(2*z.moy.inf)-1)/(exp(2*z.moy.inf)+1)
    r.com.sup <- (exp(2*z.moy.sup)-1)/(exp(2*z.moy.sup)+1)
    zeta <- 0.5*log((1+theo)/(1-theo))
    u.obs.com <- abs(z.moy-zeta)*sqrt(sum(length(var1.2),length(var3.2))-6)
    p.com <- min(pnorm(u.obs.com,0,1),pnorm(u.obs.com,0,1,lower.tail=FALSE))*2
    tab.com <- data.frame("inf"=r.com.inf,"r"=r.com,"sup"=r.com.sup,"theoretical"=theo,"U"=u.obs.com,"Pr(>|U|)"=p.com,
	" "=.psignif(p.com),stringsAsFactors=FALSE,check.names=FALSE)
    result$conf.level <- conf.level
    result$common.name <- paste("        Common correlation coefficient, ",100*conf.level,"% confidence interval\n",
	"          and equality to given value ",theo,sep="")
    result$common <- tab.com
  }
  class(result) <- if (p>alpha) {"RVtest"} else {"htest"}
  return(result)
}

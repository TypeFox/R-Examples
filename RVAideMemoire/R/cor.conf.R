cor.conf <-
function(var1,var2,theo) {
  if (length(var1)!=length(var2)) {stop("'var1' and 'var2' lengths differ")}
  if (theo<(-1) | theo>1) {stop("'theo' must be between -1 and 1")}
  nul <- as.numeric(row.names(table(c(which(is.na(var1)),which(is.na(var2))))))
  var1.2 <- if(length(nul)>0) {var1[-nul]} else {var1}
  var2.2 <- if(length(nul)>0) {var2[-nul]} else {var2}
  dname <- paste(deparse(substitute(var1))," and ",deparse(substitute(var2)),sep="")
  r <- as.numeric(cor.test(var1.2,var2.2,method="pearson")$estimate)
  names(r) <- "cor"
  z <- 0.5*log((1+r)/(1-r))
  zeta <- 0.5*log((1+theo)/(1-theo))
  u.obs <- abs(z-zeta)*sqrt(length(var1.2)-3)
  names(u.obs) <- "U"
  p <- 2*min(pnorm(u.obs,0,1),pnorm(u.obs,0,1,lower.tail=FALSE))
  met <- "Equality of a Pearson's linear correlation coefficient to a given value"
  alternative <- "two.sided"
  nval <- theo
  names(nval) <- "coefficient"
  result <- list(method=met,data.name=dname,statistic=u.obs,p.value=p,alternative=alternative,
    null.value=nval,estimate=r)
  class(result) <- "htest"
  return(result)
}

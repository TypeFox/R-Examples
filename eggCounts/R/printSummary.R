estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

printSummary<-function(m){
  summary<-t(apply(m,2,function(x) round(c(mean=mean(x), sd=sd(x), quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)),HPDLow95=HPDinterval(mcmc(cbind(reduction=x)), 0.95)[1,1],mode=estimate_mode(x),HPDHigh95=HPDinterval(mcmc(cbind(reduction=x)), 0.95)[1,2]),4)))
  rownames(summary)<-colnames(m)
  print(summary)
}
SummaryPts <-function(object, ...) UseMethod("SummaryPts")

SummaryPts.default <- function(object, mu,Sigma,alphasens = 1, alphafpr = 1,
                           n.iter = 10^6, FUN, ...){
  samples <- rmvnorm(n.iter, mu, Sigma)
  sens <- inv.trafo(alphasens,samples[,1])
  fpr <- inv.trafo(alphafpr,samples[,2])
  out <- lapply(FUN, function(x){x(sens, fpr)})
  class(out) <- "SummaryPts"
  out
}

SummaryPts.reitsma <- function(object, n.iter = 10^6, FUN = NULL, ...){
  fit <- object
  if(length(coef(fit)) > 2){
    stop("SummaryPts is not be used for meta-regression!")}
  if(is.null(FUN)){FUN <- list(posLR = function(sens,fpr){sens/fpr},
                               negLR = function(sens,fpr){(1-sens)/(1-fpr)},
                               invnegLR = function(sens, fpr){(1-fpr)/(1-sens)},
                               DOR = function(sens, fpr){sens*(1-fpr)/((1-sens)*fpr)})}
  SummaryPts.default(mu = coef(fit)[1:2], Sigma = vcov(fit), 
                         alphasens = fit$alphasens, 
                         alphafpr = fit$alphafpr,
                         n.iter = n.iter, FUN = FUN)
  }

print.SummaryPts <- function(x, ...){
print(lapply(x, mean))
}

summary.SummaryPts <- function(object, level = 0.95, digits = 3, ...){
  x <- object
  quantiles <- c((1-level)/2, 1-(1-level)/2) 
  qq <- t(sapply(x, function(x){stats::quantile(x, probs = quantiles)}))
  qq <- signif(cbind(sapply(x,mean), sapply(x,median), qq), digits)
  colnames(qq)[1:2] <- c("Mean", "Median")
  qq
}
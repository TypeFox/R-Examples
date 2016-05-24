summary.eiRegBayes <- function(object, CI = 0.95, ...){
  coef <- apply(object$draws, c(1,2), mean)
  se <- apply(object$draws, c(1,2), sd)
  cc <- apply(object$draws, c(1,2), quantile, c((1-CI)/2, 1-(1-CI)/2))
  
  quants <- matrix(cc, nrow = prod(dim(coef)[1:2]), ncol = 2, byrow = TRUE)
  nidx <- apply(expand.grid(dimnames(coef)), 1, paste, collapse = ".")
  tab <- cbind(c(coef), c(se), quants)
  colnames(tab) <- c("Mean", "Std. Dev.", rownames(cc))                     
  rownames(tab) <- nidx
  
  out <- list(call = object$call, coef = tab, sims = dim(object$draws)[3])
  class(out) <- "eiRegBayesSum"
  out
}

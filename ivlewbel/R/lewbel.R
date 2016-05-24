lewbel <-
function(formula, data, clustervar = NULL, robust = TRUE){
  res = lewbel.est(formula, data, clustervar, robust)
  est = NULL
  est$coef.est = res$coefs 
  est$call = match.call() 
  est$num.obs = res$num.obs
  est$j.test = res$j.test
  est$f.test.stats = res$f.test.stats
  class(est) = "lewbel.model"
  est  
}

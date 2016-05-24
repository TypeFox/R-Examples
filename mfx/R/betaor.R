betaor <-
function(formula, data, robust = FALSE, clustervar1 = NULL, 
                   clustervar2 = NULL, control = betareg.control(), 
                   link.phi = NULL, type = "ML"){
  
  res = betaorest(formula, data, robust, clustervar1, clustervar2, 
                  control, link.phi, type)
  
  est = NULL
  zstat = log(res$or$oddsratio)*res$or$oddsratio/res$or$se
  
  est$oddsratio = cbind(OddsRatio = res$or$oddsratio,
                        StdErr = res$or$se,
                        z.value = zstat,
                        p.value = 2*pt(-abs(zstat), df = Inf))
  
  colnames(est$oddsratio) = c("OddsRatio","Std. Err.","z","P>|z|")
  rownames(est$oddsratio) =  rownames(res$or)
  est$fit = res$fit
  est$call = match.call() 
  class(est) = "betaor"
  est
}

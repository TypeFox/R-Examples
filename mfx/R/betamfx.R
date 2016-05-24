betamfx <-
function(formula, data, atmean = TRUE, robust = FALSE, clustervar1 = NULL, 
                   clustervar2 = NULL, control = betareg.control(), 
                   link.phi = NULL, type = "ML"){
  
  res = betamfxest(formula = formula, data = data, atmean, robust, 
                   clustervar1 = clustervar1, clustervar2 = clustervar2,
                   control = control, link.phi = link.phi, type = type)
  
  est = NULL
  est$mfxest = cbind(dFdx = res$mfx$mfx,
                     StdErr = res$mfx$se,
                     z.value = res$mfx$mfx/res$mfx$se,
                     p.value = 2*pt(-abs(res$mfx$mfx/res$mfx$se), df = Inf))
  colnames(est$mfxest) = c("dF/dx","Std. Err.","z","P>|z|")
  rownames(est$mfxest) =  rownames(res$mfx)
  
  est$fit = res$fit
  est$dcvar = rownames(res$mfx[res$mfx$discretechgvar==1,])  
  est$call = match.call() 
  class(est) = "betamfx"
  est
}

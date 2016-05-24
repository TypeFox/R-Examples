negbinmfx <-
function(formula, data, atmean=TRUE, robust=FALSE, clustervar1=NULL, 
                     clustervar2=NULL, start = NULL, control = glm.control()){
  
  res = negbinmfxest(formula, data, atmean, robust, clustervar1, 
                     clustervar2, start = start, control = control)
  
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
  class(est) = "negbinmfx"
  est
}

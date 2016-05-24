negbinirr <-
function(formula, data, robust = FALSE, clustervar1=NULL, 
                     clustervar2=NULL, start = NULL, control = glm.control()){

  res = negbinirrest(formula, data, robust, clustervar1, clustervar2, 
                     start = start, control = control)
  
  est = NULL
  zstat = log(res$irr$irr)*res$irr$irr/res$irr$se
  
  est$irr = cbind(IRR = res$irr$irr,
                        StdErr = res$irr$se,
                        z.value = zstat,
                        p.value = 2*pt(-abs(zstat), df = Inf))
  colnames(est$irr) = c("IRR","Std. Err.","z","P>|z|")
  rownames(est$irr) =  rownames(res$irr)
  
  est$fit = res$fit
  est$call = match.call() 
  class(est) = "negbinirr"
  est
}

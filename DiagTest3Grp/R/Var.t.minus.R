Var.t.minus <-
function(mu.minus,mu0,s.minus,s0,n.minus,n0)
  {
    ###This function calculates the variance of the estimate on the lower optimal cut-point,under normality assumption
    ###To get the counterparts for t.plus on mu0,mu.plus,s0,s.plus, simply simultaneously replace mu.minus by mu0, mu0 by mu.plus, s.minus by s0 and s0 by s.plus
    
    res1 <- PartialDeriv.optimalCutoff(mu.minus,mu0,s.minus,s0)
    deriv.tminus.mu.minus <- res1$deriv.mu1
    deriv.tminus.mu0 <- res1$deriv.mu2
    
    deriv.tminus.s.minus <- res1$deriv.s1
    deriv.tminus.s0 <- res1$deriv.s2

    var0 <- deriv.tminus.mu.minus^2*var.mu(s.minus,n.minus)+deriv.tminus.mu0^2*var.mu(s0,n0)+deriv.tminus.s.minus^2*var.sigma(s.minus,n.minus)+deriv.tminus.s0^2*var.sigma(s0,n0)
    return(var0)    
  }


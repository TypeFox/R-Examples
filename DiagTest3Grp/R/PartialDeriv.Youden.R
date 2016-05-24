PartialDeriv.Youden <-
function(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus)
  {
    #####This function calculates the partial derivatives of the Youden index w.r.t the normal parameters (the 3 means and 3 SDs) and the two cut-point
    ###t.minus partial derivatives w.r.t relevant parameters mu.minus and mu0, s.minus, s0
    res1 <- PartialDeriv.optimalCutoff(mu.minus,mu0,s.minus,s0)
    deriv.tminus.mu.minus <- res1$deriv.mu1
    deriv.tminus.mu0 <- res1$deriv.mu2
    
    deriv.tminus.s.minus <- res1$deriv.s1
    deriv.tminus.s0 <- res1$deriv.s2

    ###t.plus partial derivatives w.r.t relevant parameters:mu0,mu.plus,s0,s.plus
    res2 <- PartialDeriv.optimalCutoff(mu0,mu.plus,s0,s.plus)
    deriv.tplus.mu0 <- res2$deriv.mu1
    deriv.tplus.mu.plus <- res2$deriv.mu2
    
    deriv.tplus.s0 <- res2$deriv.s1
    deriv.tplus.s.plus <- res2$deriv.s2
    
    ####partial derivatives of Youden index Y on all normal parameters (note that first Y=Se+Sp+Sm-1, omitting the 0.5 multiplication factor here which will be multiplied at the end)
    stand.tminus <- (t.minus-mu.minus)/s.minus

    #Y.mu.minus <- deriv.tminus.mu.minus*(1/s.minus*dnorm(stand.tminus)-1/s0*dnorm((t.minus-mu0)/s0))-dnorm(stand.tminus)
    Y.mu.minus <- deriv.tminus.mu.minus*(1/s.minus*dnorm(stand.tminus)-1/s0*dnorm((t.minus-mu0)/s0))-1/s.minus*dnorm(stand.tminus)
    
    stand.tplus <- (t.plus-mu.plus)/s.plus
    
    #Y.mu.plus <- deriv.tplus.mu.plus*(1/s0*dnorm((t.plus-mu0)/s0)-1/s.plus*dnorm(stand.tplus))+dnorm(stand.tplus)
    Y.mu.plus <- deriv.tplus.mu.plus*(1/s0*dnorm((t.plus-mu0)/s0)-1/s.plus*dnorm(stand.tplus))+1/s.plus*dnorm(stand.tplus)
    
    #Y.mu0 <- deriv.tminus.mu0*( 1/s.minus*dnorm(stand.tminus)-1/s0*dnorm((t.minus-mu0)/s0))+deriv.tplus.mu0*(1/s0*dnorm((t.plus-mu0)/s0)-1/s.plus*dnorm(stand.tplus))+dnorm((t.minus-mu0)/s0)-dnorm((t.plus-mu0)/s0)
    Y.mu0 <- deriv.tminus.mu0*( 1/s.minus*dnorm(stand.tminus)-1/s0*dnorm((t.minus-mu0)/s0))+deriv.tplus.mu0*(1/s0*dnorm((t.plus-mu0)/s0)-1/s.plus*dnorm(stand.tplus))+1/s0*(dnorm((t.minus-mu0)/s0)-dnorm((t.plus-mu0)/s0))
    
    Y.s.minus <- deriv.tminus.s.minus*( 1/s.minus*dnorm(stand.tminus)-1/s0*dnorm((t.minus-mu0)/s0))-stand.tminus/s.minus*dnorm(stand.tminus)
    
    Y.s.plus <- deriv.tplus.s.plus*(1/s0*dnorm((t.plus-mu0)/s0)-1/s.plus*dnorm(stand.tplus))+stand.tplus/s.plus*dnorm(stand.tplus)
    
    Y.s0 <- deriv.tminus.s0*( 1/s.minus*dnorm(stand.tminus)-1/s0*dnorm((t.minus-mu0)/s0))+deriv.tplus.s0*(1/s0*dnorm((t.plus-mu0)/s0)-1/s.plus*dnorm(stand.tplus))+(t.minus-mu0)/(s0^2)*dnorm((t.minus-mu0)/s0)-(t.plus-mu0)/(s0^2)*dnorm((t.plus-mu0)/s0)

    ###Y=0.5*(Se+Sp+Sm-1), the derivatives omit the 0.5 factor
    Y.mu.minus <- 0.5*Y.mu.minus
    Y.mu.plus <- 0.5*Y.mu.plus
    Y.mu0 <- 0.5*Y.mu0

    Y.s.minus <- 0.5*Y.s.minus
    Y.s.plus <- 0.5*Y.s.plus
    Y.s0 <- 0.5*Y.s0

    return(list(Y.mu.minus=Y.mu.minus,Y.mu.plus=Y.mu.plus,Y.mu0=Y.mu0,Y.s.minus=Y.s.minus,Y.s.plus=Y.s.plus,Y.s0=Y.s0))
  }


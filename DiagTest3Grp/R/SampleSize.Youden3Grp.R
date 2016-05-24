SampleSize.Youden3Grp <-
function(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus,lam.minus=1/3,lam0=1/3,lam.plus=1/3,typeIerror=0.05,margin=0.05)
  {
    #####This function calculates the sample size to collect data on one biomarker in order to estimate the Youden index of the biomarker within margin of error
    ##########Input:
    ###mu.minus,mu0,mu.plus: normal means, can be specified by users arbitrarily (must be in increasing order) or estimated from data first
    ###s.minus,s0,s.plus:normal standard deviations,can be specified by users or estimated from data first
    ###t.minus,t.plus: optimal cut points
    ###lam.minus, lam0,lam.plus
    ###FisherZ: For sample size calculation purpose, place margin of error on the CI of the Fisher's Z transformed Youden estimate
    ###youden: default to NULL but can be provided if estimation of youden is known, for sample size calculation under Fisher's Z transformation
    J.partial <- PartialDeriv.Youden(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus)
    Mj <- s.minus^2*lam0/lam.minus*(J.partial$Y.mu.minus^2+0.5*J.partial$Y.s.minus^2)+s0^2*(J.partial$Y.mu0^2+0.5*J.partial$Y.s0^2)+s.plus^2*lam0/lam.plus*(J.partial$Y.mu.plus^2+0.5*J.partial$Y.s.plus^2)
    z0 <- qnorm(typeIerror/2,lower.tail=F)    

    
    sampleSize <- z0^2*Mj/margin^2
          
    #sampleSize <- ceiling(sampleSize)
    
    return(sampleSize)
  }


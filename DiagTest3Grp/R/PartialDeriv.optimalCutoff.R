PartialDeriv.optimalCutoff <-
function(mu.minus,mu0,s.minus,s0)
  {
    ###This function calculates under normal method the partial derivatives of the lower optimal cut-point t.minus w.r.t the normal parameters (mu.minus,mu0,s.minus,s0), see Reference paper Luo&Xiong 2011
    ##To obtain the counterparts for t.plus on mu0,mu.plus,s0,s.plus, simply simultaneously replace mu.minus by mu0, mu0 by mu.plus, s.minus by s0 and s0 by s.plus

    if(s.minus==s0)
      {
        deriv.mu.minus <- 0.5
        deriv.mu0 <- 0.5
        deriv.s.minus <- NA
        deriv.s0 <- NA
      }        
    else
      {
        var.ratio <- 2*(log(s.minus)-log(s0))
        b <- s.minus^2-s0^2
        mu.diff <-mu.minus-mu0
        
        delta <- mu.diff^2+b*var.ratio##the term under sqr root in the expression of t.minus
        c <- mu0*s.minus^2-mu.minus*s0^2-s.minus*s0*sqrt(delta)###the numerator in the expression of t.minus
        
        a1 <- 2*mu0*s.minus-s0*sqrt(delta)-s0/sqrt(delta)*(s.minus^2*var.ratio+b)##derivative of the numerator in t_ on s.minus
        a2 <- -2*mu.minus*s0-s.minus*sqrt(delta)+s.minus/sqrt(delta)*(s0^2*var.ratio+b)##derivative of the numerator in t_ on s0
        
        deriv.mu.minus <- -s0/b*(s0+s.minus/sqrt(delta)*mu.diff)
        deriv.mu0 <- s.minus/b*(s.minus+s0/sqrt(delta)*mu.diff)
        deriv.s.minus <- (a1*b-2*s.minus*c)/(b^2)        
        deriv.s0 <- (a2*b+2*s0*c)/(b^2)
      }
    
    return(data.frame(deriv.mu1=deriv.mu.minus,deriv.mu2=deriv.mu0,deriv.s1=deriv.s.minus,deriv.s2=deriv.s0))
  }


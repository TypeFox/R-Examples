Var.Youden <-
function(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,n.minus,n0,n.plus,t.minus,t.plus)
  {
    ###This function calculates the Youden index's associated variance under normality , called by Youden3Grp.Variance.Normal.R


    res0 <- PartialDeriv.Youden(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus)
    Y.mu.minus <- res0$Y.mu.minus##derivative of youden w.r.t mu.minus
    Y.mu.plus <- res0$Y.mu.plus
    Y.mu0 <- res0$Y.mu0
    Y.s.minus <- res0$Y.s.minus
    Y.s.plus <- res0$Y.s.plus
    Y.s0 <- res0$Y.s0
    
    var0 <- Y.mu.minus^2*var.mu(s.minus,n.minus)+Y.mu.plus^2*var.mu(s.plus,n.plus)+Y.mu0^2*var.mu(s0,n0)+Y.s.minus^2*var.sigma(s.minus,n.minus)+Y.s.plus^2*var.sigma(s.plus,n.plus)+Y.s0^2*var.sigma(s0,n0)
    
    #print(c(n.minus=n.minus,n0=n0,n.plus=n.plus,s.minus=s.minus,s0=s0,s.plus=s.plus,Y.mu.minus=Y.mu.minus,Y.mu0=Y.mu0,Y.mu.plus=Y.mu.plus,Y.s.minus=Y.s.minus,Y.s0=Y.s0,Y.s.plus=Y.s.plus))
    #write.table(data.frame(n.minus=n.minus,n0=n0,n.plus=n.plus,mu.minus=mu.minus,mu0=mu0,mu.plus=mu.plus,Y.mu.minus=Y.mu.minus,Y.mu0=Y.mu0,Y.mu.plus=Y.mu.plus,s.minus=s.minus,s0=s0,s.plus=s.plus,Y.mu.minus=Y.mu.minus,Y.mu0=Y.mu0,Y.mu.plus=Y.mu.plus,Y.s.minus=Y.s.minus,Y.s0=Y.s0,Y.s.plus=Y.s.plus),file="EachMarker.YoudenVarianceComponents.csv",sep=",",col.names=T,row.names=F,append=T)
    return(list(var0=var0,partialDeriv=data.frame(s.minus=s.minus,s0=s0,s.plus=s.plus,Y.mu.minus=Y.mu.minus,Y.mu0=Y.mu0,Y.mu.plus=Y.mu.plus,Y.s.minus=Y.s.minus,Y.s0=Y.s0,Y.s.plus=Y.s.plus)))
  }


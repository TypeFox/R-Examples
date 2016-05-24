expandbound<-function(g, pct=.10)
  {
    ###  expand the bounds of a vector range by a specified percent
    r = range(g)
    dg = pct*diff(r)
    
    return(c(r[1]-dg, r[2]+dg)) 
  }

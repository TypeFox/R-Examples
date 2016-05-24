Phinv <-
  function(p,l,u){
    # computes with precision the quantile function
    # of the standard normal distribution,
    # truncated to the interval [l,u], using erfcinv.
    I=u<0;l[I]=-l[I];u[I]=-u[I]; # use symmetry of normal
    pu=pnorm(u,lower.tail = FALSE);
    pl=pnorm(l,lower.tail = FALSE);
    x=qnorm(pl+(pu-pl)*p,lower.tail = FALSE);
    x[I]=-x[I]; # adjust sign due to symmetry
    return(x)
  }

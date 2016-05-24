tn <-
  function(l,u)
  { # samples a column vector of length=length(l)=length(u)
    # from the standard multivariate normal distribution,
    # truncated over the region [l,u], where -a<l<u<a for some
    # 'a' and l and u are column vectors;
    # uses acceptance rejection and inverse-transform method;
    tol=2.05 # controls switch between methods
    # threshold can be tuned for maximum speed for each platform
    # case: abs(u-l)>tol, uses accept-reject from randn
    x=l;I=(abs(u-l)>tol)
    if (any(I)){
      tl=l[I];tu=u[I];x[I]=trnd(tl,tu)
    }
    # case: abs(u-l)<tol, uses inverse-transform
    I<-!I
    if (any(I)){
      tl=l[I];tu=u[I];pl=pnorm(tl);pu=pnorm(tu);
      x[I]=qnorm(pl+(pu-pl)*runif(length(tl)))
    }
    return(x)
  }

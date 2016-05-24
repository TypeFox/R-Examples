newton <-
  function(p,l,u){
    ql=qfun(l);
    qu=rep(0,length(u)); # deal with infinite limits
    I=is.infinite(u);
    I=!I; # find where it is finite
    if (any(I)){
      qu[I]=qfun(u[I]);
    }
    l=l^2;u=u^2;
    # set initial value for Newton iteration
    x=sqrt(l-2*log(1+p*expm1(l/2-u/2)))
    # initialize Newton method
    err=Inf;
    while (err>10^-10){
      del=-qfun(x)+(1-p)*exp(.5*(x^2-l))*ql+p*exp(.5*(x^2-u))*qu
      x=x-del # Newton's step
      err=max(abs(del)) # find the maximum error
    }
    return(x)
  }

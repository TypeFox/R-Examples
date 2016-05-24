trandn <-
  function(l,u){
    ## truncated normal generator
    # * efficient generator of a vector of length(l)=length(u)
    # from the standard multivariate normal distribution,
    # truncated over the region [l,u];
    # infinite values for 'u' and 'l' are accepted;
    # * Remark:
    # If you wish to simulate a random variable
    # 'Z' from the non-standard Gaussian N(m,s^2)
    # conditional on l<Z<u, then first simulate
    # X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;
    #
    # Reference:
    # Z. I. Botev (2015),
    # "The Normal Law Under Linear Restrictions:
    #  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
    if (any(l>u)){
      stop('Truncation limits have to be vectors of the same length with l<u')
    }
    x=rep(0,length(l));
    a=.4; # treshold for switching between methods
    # threshold can be tuned for maximum speed for each Matlab version
    # three cases to consider:
    # case 1: a<l<u
    I=l>a;
    if (any(I)){
      tl=l[I]; tu=u[I]; x[I]=ntail(tl,tu);
    }
    # case 2: l<u<-a
    J=u<(-a);
    if (any(J)){
      tl=-u[J]; tu=-l[J]; x[J]=-ntail(tl,tu);
    }
    # case 3: otherwise use inverse transform or accept-reject
    I=!(I|J);
    if (any(I)){
      tl=l[I]; tu=u[I]; x[I]=tn(tl,tu);
    }
    return(x)
  }

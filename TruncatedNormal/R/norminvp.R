norminvp <-
  function(p,l,u){
    ## normal quantile function with precision
    # computes with tail-precision the quantile function
    # of the standard normal distribution at 0<=p<=1,
    # and truncated to the interval [l,u];
    # Inf values for vectors 'l' and 'u' accepted;
    #
    # * Example 1:
    # # Suppose you wish to simulate a random variable
    # # 'Z' from the non-standard Gaussian N(m,s^2)
    # # conditional on l<Z<u. First compute
    #  m=1;l=10;u=20;s=1; 
    #  X=norminvp(runif(1),(l-m)/s,(u-m)/s); # and then set
    #  Z=m+s*X
    #
    # Reference:
    # Z. I. Botev (2015),
    # "The Normal Law Under Linear Restrictions:
    #  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
    
    
    if ((length(l)!=length(p))|any(l>u)|any(p>1)|any(p<0)){
      stop('l, u, and p must be the same length with u>l and 0<=p<=1')
    }
    x=rep(NaN,length(l)); # allocate memory
    I=(p==1);x[I]=u[I]; # extreme values of quantile
    J=(p==0);x[J]=l[J];
    I=!(I|J); # cases for which 0<x<1
    if (any(I)){
      x[I]=cases(p[I],l[I],u[I]);
    }
    return(x)
  }

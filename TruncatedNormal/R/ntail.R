ntail <-
  function(l,u)
  {# samples a column vector of length=length(l)=length(u)
    # from the standard multivariate normal distribution,
    # truncated over the region [l,u], where l>0 and
    # l and u are column vectors;
    # uses acceptance-rejection from Rayleigh distr;
    # method due to Marsaglia (1964);
    c=l^2/2;n=length(l);f=expm1(c-u^2/2)
    x=c-log(1+runif(n)*f); # sample using Rayleigh
    # keep list of rejected
    I=which(runif(n)^2*x>c);d=length(I)
    while (d>0) { # while there are rejections
      cy=c[I] # find the thresholds of rejected
      y=cy-log(1+runif(d)*f[I])
      idx=runif(d)^2*y<cy # accepted
      x[I[idx]]=y[idx] # store the accepted
      I=I[!idx] # remove accepted from list
      d=length(I) # number of rejected
    }
    x=sqrt(2*x) # this Rayleigh transform can be delayed till the end
  }

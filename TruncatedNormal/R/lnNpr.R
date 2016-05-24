lnNpr <-
  function(a,b)
  { # computes ln(P(a<Z<b))
    # where Z~N(0,1) very accurately for any 'a', 'b'
    p=rep(0,length(a))
    # case b>a>0
    I=a>0
    if (any(I)){
      pa=pnorm(a[I],lower.tail = FALSE, log.p = TRUE)
      pb=pnorm(b[I],lower.tail = FALSE, log.p = TRUE)
      p[I]=pa+log1p(-exp(pb-pa)) 
    }
    # case a<b<0
    idx=b<0
    if (any(idx)){
      pa=pnorm(a[idx], log.p = TRUE)
      pb=pnorm(b[idx], log.p = TRUE)
      p[idx]=pb+log1p(-exp(pa-pb))
    }
    # case a<0<b
    I=!I&!idx
    if (any(I)){
      pa=pnorm(a[I])
      pb=pnorm(b[I],lower.tail = FALSE)
      p[I]=log1p(-pa-pb)
    }
    return(p)
  }

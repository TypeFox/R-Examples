mvnpr <-
  function(n,L,l,u,mu){
    # computes P(l<X<u), where X is normal with
    # 'Cov(X)=L*L' and zero mean vector;
    # exponential tilting uses parameter 'mu';
    # Monte Carlo uses 'n' samples;
    d=length(l); # Initialization
    mu[d]=0;
    Z=matrix(0,d,n); # create array for variables
    p=0;
    for (k in 1:(d-1)){
      # compute matrix multiplication L*Z
      col=t(L[k,1:k])%*%Z[1:k,];
      # compute limits of truncation
      tl=l[k]-mu[k]-col;
      tu=u[k]-mu[k]-col;
      #simulate N(mu,1) conditional on [tl,tu]
      Z[k,]=mu[k]+trandn(tl,tu);
      # update likelihood ratio
      p = p+lnNpr(tl,tu)+.5*mu[k]^2-mu[k]*Z[k,];
    }
    # deal with final Z(d) which need not be simulated
    col=L[d,]%*%Z;tl=l[d]-col;tu=u[d]-col;
    p=p+lnNpr(tl,tu); # update LR corresponding to Z(d)
    p=exp(p); # now switch back from logarithmic scale
    prob=mean(p);relErr=sd(p)/sqrt(n)/prob; # relative error
    est=list(prob=prob,relErr=relErr)
  }

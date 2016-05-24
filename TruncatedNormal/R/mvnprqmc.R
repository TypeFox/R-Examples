mvnprqmc <-
  function(n,L,l,u,mu){
    # computes P(l<X<u), where X is normal with
    # 'Cov(X)=L*L' and zero mean vector;
    # exponential tilting uses parameter 'mu';
    # Quasi Monte Carlo uses 'n' samples;
    d=length(l); # Initialization
    mu[d]=0;
    Z=matrix(0,d,n); # create array for variables
    # QMC pointset
    x=sobol(n, dim = d-1, init =TRUE, scrambling = 1, seed=ceiling(10^6*runif(1)))
    p=0;
    for (k in 1:(d-1)){
      # compute matrix multiplication L*Z
      col=t(L[k,1:k])%*%Z[1:k,];
      # compute limits of truncation
      tl=l[k]-mu[k]-col;
      tu=u[k]-mu[k]-col;
      #simulate N(mu,1) conditional on [tl,tu] via QMC
      if (d>2){
        Z[k,]=mu[k]+norminvp(x[,k],tl,tu);
      } else {
        Z[k,]=mu[k]+norminvp(x,tl,tu);
      }
      # update likelihood ratio
      p = p+lnNpr(tl,tu)+.5*mu[k]^2-mu[k]*Z[k,];
    }
    # deal with final Z(d) which need not be simulated
    col=L[d,]%*%Z;tl=l[d]-col;tu=u[d]-col;
    p=p+lnNpr(tl,tu); # update LR corresponding to Z(d)
    p=mean(exp(p)); # now switch back from logarithmic scale
  }

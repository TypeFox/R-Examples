mvnrnd <-
  function(n,L,l,u,mu){
    # generates the proposals from the exponentially tilted 
    # sequential importance sampling pdf;
    # output:    'logpr', log-likelihood of sample
    #             Z, random sample 
    d=length(l); # Initialization
    mu[d]=0;
    Z=matrix(0,d,n); # create array for variables
    p=0;
    for (k in 1:d){
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
    return(list(logpr=p,Z=Z))
  }

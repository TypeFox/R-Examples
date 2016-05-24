cholperm <-
  function( Sig, l, u ){
    #  Computes permuted lower Cholesky factor L for Sig
    #  by permuting integration limit vectors l and u.
    #  Outputs perm, such that Sig(perm,perm)=L%*%t(L).
    #
    # Reference: 
    #  Gibson G. J., Glasbey C. A., Elston D. A. (1994),
    #  "Monte Carlo evaluation of multivariate normal integrals and
    #  sensitivity to variate ordering", 
    #  In: Advances in Numerical Methods and Applications, pages 120--126
    eps=10^-10; # round-off error tolerance
    d=length(l);perm=1:d; # keep track of permutation
    L=matrix(0,d,d);z=rep(0,d);
    for (j in 1:d){
      pr=rep(Inf,d); # compute marginal prob.
      I=j:d; # search remaining dimensions
      D=diag(Sig);
      if (j >2){
        s=D[I]-L[I,1:(j-1)]^2%*%rep(1,j-1)
      } else if (j==2){
        s=D[I]-L[I,1]^2
      } else {
        s=D[I]
      }
      s[s<0]=eps;s=sqrt(s);
      if (j >2){
        cols=L[I,1:(j-1)]%*%z[1:(j-1)]
      } else if (j==2){
        cols=L[I,1]*z[1]
      } else {
        cols=0
      }
      tl=(l[I]-cols)/s;tu=(u[I]-cols)/s;pr[I]=lnNpr(tl,tu);
      # find smallest marginal dimension
      k=which.min(pr);
      # flip dimensions k-->j
      jk=c(j,k);kj=c(k,j);
      Sig[jk,]=Sig[kj,];Sig[,jk]=Sig[,kj]; # update rows and cols of Sig
      L[jk,]=L[kj,]; # update only rows of L
      l[jk]=l[kj];u[jk]=u[kj]; # update integration limits
      perm[jk]=perm[kj]; # keep track of permutation
      # construct L sequentially via Cholesky computation
      s=Sig[j,j]-sum(L[j,1:(j-1)]^2);
      if (s<(-0.001)){
        stop('Sigma is not positive semi-definite')
      }
      s[s<0]=eps;L[j,j]=sqrt(s);
      if (j<d){
        if (j >2){
          L[(j+1):d,j]=(Sig[(j+1):d,j]-L[(j+1):d,1:(j-1)]%*%L[j,1:(j-1)])/L[j,j];
        } else if (j==2){
          L[(j+1):d,j]=(Sig[(j+1):d,j]-L[(j+1):d,1]*L[j,1])/L[j,j];
        } else if (j==1){
          L[(j+1):d,j]=Sig[(j+1):d,j]/L[j,j];
        }
      }
      # find mean value, z(j), of truncated normal:
      tl=(l[j]-L[j,1:j]%*%z[1:j])/L[j,j];
      tu=(u[j]-L[j,1:j]%*%z[1:j])/L[j,j];
      w=lnNpr(tl,tu); # aids in computing expected value of trunc. normal
      z[j]=(exp(-.5*tl^2-w)-exp(-.5*tu^2-w))/sqrt(2*pi);
    }
    return(list(L=L,l=l,u=u,perm=perm))
  }

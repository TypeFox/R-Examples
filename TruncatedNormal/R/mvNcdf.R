mvNcdf <-
  function(l,u,Sig,n){
    ## truncated multivariate normal cumulative distribution
    # computes an estimator of the probability Pr(l<X<u),
    # where 'X' is a zero-mean multivariate normal vector
    # with covariance matrix 'Sig', that is, X~N(0,Sig)
    # infinite values for vectors 'u' and 'l' are accepted;
    # Monte Carlo method uses sample size 'n'; the larger
    # the 'n', the smaller the relative error of the estimator;
    #      output: structure 'est' with
    #              1. estimated value of probability Pr(l<X<u)
    #              2. estimated relative error of estimator
    #              3. theoretical upper bound on true Pr(l<X<u)
    #              Remark: If you want to estimate Pr(l<Y<u),
    #                   where Y~N(m,Sig) has mean vector 'm',
    #                     then use 'mvNcdf(Sig,l-m,u-m,n)'.
    # * Example:
    #  d=25;l=rep(5,d);u=rep(Inf,d);
    #  Sig=0.5*diag(d)+.5*matrix(1,d,d);
    #  est=mvNcdf(l,u,Sig,10^4) # output of our method
    #
    # Reference:
    # Z. I. Botev (2015),
    # "The Normal Law Under Linear Restrictions:
    #  Simulation and Estimation via Minimax Tilting",
    #  submitted to JRSS(B)
    d=length(l); # basic input check
    if  (length(u)!=d|d!=sqrt(length(Sig))|any(l>u)){
      stop('l, u, and Sig have to match in dimension with u>l')
    }
    # Cholesky decomposition of matrix
    out=cholperm( Sig, l, u ); L=out$L; l=out$l; u=out$u; D=diag(L);
    if (any(D<10^-10)){
      warning('Method may fail as covariance matrix is singular!')
    }
    L=L/D;u=u/D;l=l/D; # rescale
    L=L-diag(d); # remove diagonal
    # find optimal tilting parameter via non-linear equation solver
    xmu<-nleq(l,u,L) # nonlinear equation solver 
    #soln<-nleqslv(rnorm(2*d-2),gradpsi,jac=NULL,L,l,u,method = c("Broyden"))
    #xmu<-soln$x;
    x=xmu[1:(d-1)];mu=xmu[d:(2*d-2)]; # assign saddlepoint x* and mu*
    est=mvnpr(n,L,l,u,mu);
    # compute psi star
    est$upbnd=exp(psy(x,L,l,u,mu));
    return(est)
  }

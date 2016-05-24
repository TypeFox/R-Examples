mvNqmc <-function(l,u,Sig,n){
  ## truncated multivariate normal cumulative distribution (qmc version)
  # computes an estimator of the probability Pr(l<X<u),
  # where 'X' is a zero-mean multivariate normal vector
  # with covariance matrix 'Sig', that is, X~N(0,Sig)
  # infinite values for vectors 'u' and 'l' are accepted;
  # 
  # This version uses a Quasi Monte Carlo (QMC) pointset
  # of size ceil(n/12) and estimates the relative error
  # using 12 independent randomized QMC estimators; QMC
  # is slower than ordinary Monte Carlo (see my mvncdf.m), 
  # but is also likely to be more accurate when d<50.
  #
  # output:      structure 'est' with
  #              1. estimated value of probability Pr(l<X<u)
  #              2. estimated relative error of estimator
  #              3. theoretical upper bound on true Pr(l<X<u)
  #              
  # * Remark: If you want to estimate Pr(l<Y<u),
  #           where Y~N(m,Sig) has mean vector 'm',
  #           then use 'mvNqmc(Sig,l-m,u-m,n)'.
  #
  # * Example:
  #  d=25;l=rep(5,d);u=rep(Inf,d);
  #  Sig=0.5*diag(d)+.5*matrix(1,d,d);
  #  est=mvNqmc(l,u,Sig,10^4) # output of our method
  #
  # Reference: Z. I. Botev (2015),
  # "The Normal Law Under Linear Restrictions:
  #  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
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
  x=xmu[1:(d-1)];mu=xmu[d:(2*d-2)]; # assign saddlepoint x* and mu*
  p=rep(NaN,12);
  for (i in 1:12){ # repeat randomized QMC
    p[i]=mvnprqmc(ceiling(n/12),L,l,u,mu);
  }
  prob=mean(p); # average of QMC estimates
  relErr=sd(p)/sqrt(12)/prob; # relative error
  upbnd=exp(psy(x,L,l,u,mu)); # compute psi star
  est=list(prob=prob,relErr=relErr,upbnd=upbnd)
}

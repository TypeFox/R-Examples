mvrandn <-
  function(l,u,Sig,n){
    ## truncated multivariate normal generator
    # simulates 'n' random vectors exactly/perfectly distributed
    # from the d-dimensional N(0,Sig) distribution (zero-mean normal
    # with covariance 'Sig') conditional on l<X<u;
    # infinite values for 'l' and 'u' are accepted;
    # output:   'd' times 'n' array 'rv' storing random vectors;
    #
    # * Example:
    #  d=60;n=10^3;Sig=0.9*matrix(1,d,d)+0.1*diag(d);l=(1:d)/d*4;u=l+2;
    #  X=mvrandn(l,u,Sig,n);boxplot(t(X)) # plot marginals
    #
    # * Notes: Algorithm may not work if 'Sig' is close to being rank deficient.
    # Reference:
    # Z. I. Botev (2015), "The Normal Law Under Linear Restrictions:
    #  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
    
    d=length(l); # basic input check
    if  (length(u)!=d|d!=sqrt(length(Sig))|any(l>u)){
      stop('l, u, and Sig have to match in dimension with u>l')
    }
    # Cholesky decomposition of matrix
    out=cholperm(Sig,l,u);
    Lfull=out$L;l=out$l;u=out$u;D=diag(Lfull);perm=out$perm;
    if (any(D<10^-10)){
      warning('Method may fail as covariance matrix is singular!')
    }
    L=Lfull/D;u=u/D;l=l/D; # rescale
    L=L-diag(d); # remove diagonal
    # find optimal tilting parameter via non-linear equation solver
    xmu<-nleq(l,u,L) # nonlinear equation solver 
    x=xmu[1:(d-1)];mu=xmu[d:(2*d-2)]; # assign saddlepoint x* and mu*
    # compute psi star
    psistar=psy(x,L,l,u,mu);
    # start acceptance rejection sampling
    iter=0; rv=c();
    repeat{ 
      out=mvnrnd(n,L,l,u,mu);logpr=out$logpr;Z=out$Z; # simulate n proposals
      idx=-log(runif(n))>(psistar-logpr); # acceptance tests
      rv=cbind(rv,Z[,idx]);  # accumulate accepted
      accept=dim(rv)[2]; # keep track of # of accepted
      iter=iter+1;  # keep track of while loop iterations
      if (iter==10^3){ # if iterations are getting large, give warning
        warning('Acceptance prob. smaller than 0.001')
      } else if (iter>10^4){ # if iterations too large, seek approximation only
        accept=n;rv=cbind(rv,Z); # add the approximate samples
        warning('Sample is only approximately distributed.')
      }
      if (accept>=n){# if # of accepted is less than n
        break
      }
    }
    # finish sampling; postprocessing
    out=sort(perm,decreasing = FALSE,index.return = TRUE);order=out$ix;
    rv=rv[,1:n]; # cut-down the array to desired n samples
    rv=Lfull%*%rv; # reverse scaling of L
    rv=rv[order,]; # reverse the Cholesky permutation
  }

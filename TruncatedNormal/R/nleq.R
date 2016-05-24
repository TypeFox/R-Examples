nleq <-
  function(l,u,L){
    d=length(l);x=rep(0,2*d-2) # initial point for Newton iteration
    err=Inf; iter=0;
    while (err>10^-10){
      f=gradpsi(x,L,l,u)
      Jac=f$Jac
      grad=f$grad
      del=solve(Jac,-grad) # Newton correction
      x=x+del
      err=sum(grad^2)
      iter=iter+1
      if (iter>100){
        stop('Covariance matrix is ill-conditioned and method failed')
      }
    }
    return(x)
  }

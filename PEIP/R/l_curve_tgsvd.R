l_curve_tgsvd <-
function(U,d,X,Lam,G,L)
  {
### %       [rho,eta,reg_param,m] = l_curve_tgsvd
### % Initialization.
    di=dim(G);

    M = di[1]
    N = di[2]

    p=rnk(L);

    reg_param=  (1:(p+1))

    Y= t( solve(X) ) ;

    lambda=sqrt(diag( Lam)  * diag(Lam ));

    dproj=rep(0, M)
    
###for(i in 1:M)
###  {
###dproj[i] = t(U[,i]) %*% d;
###}

    dproj= t(U)  %*% d;
                 

    m=matrix( rep(0, length=N*M), nrow=N, ncol=M)
             
### %build the solutions
    iend = length(d)
    
    m[,1]=(dproj[iend]/lambda[iend]) * Y[,iend]
             
    for(q in 2:N)
      {
        m[,q]=m[,q-1]+(dproj[iend-q+1]/lambda[iend-q+1]) * Y[,iend-q+1];
      }

### % initialize output variables
    eta = rep(0, M) ;
    rho = eta; 

### % calculate the solution misfit and seminorm for each generalized singular value used
    for(i in 1:M)
      {
        rho[i] =  Mnorm(G %*% m[,i]-d);
        eta[i] =  Mnorm(L %*% m[,i]);
      }

    return(list( rho=rho,eta=eta,reg_param=reg_param,m=m    ))
    

  }

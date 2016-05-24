l_curve_tikh_gsvd <-
function(U,d,X,Lam,Mu,G,L,npoints, varargin=NULL )
  {
##########  [rho,eta,reg_param,m] = l_curve_tikh_gsvd
    eps = .Machine$double.eps

### % Initialization.
    di=dim(G);

    M = di[1]
    N = di[2]

    p=rnk(L);

    Y= t( solve(X) ) ;


    lambda=sqrt(diag( Lam)  * diag(Lam ));


    mu=sqrt(diag(Mu)*diag(Mu));


    dproj= U %*% d;

    MYmu = c(mu, rep(0, times=length(lambda)-length(mu)) )
    gamma=lambda/ MYmu;

    m=  matrix(rep(0, N*npoints), ncol=npoints)

    rho= rep(0, npoints);
    eta=rep(0, npoints );

###

    if(length(varargin)==0)
      {
        smin_ratio = 16*eps;

### % Minimum regularization parameter
        reg_param = rep(0, length=npoints)
        
        if(M<=N)
          {
            reg_param[npoints] = max(c(gamma[p],gamma[N-M+1]*smin_ratio)  ); 

### % ratio so that reg_param(1) will be s(1)
            ratio = (gamma[N-M+1]/reg_param[npoints])^(1/(npoints-1));
            
          }
        else
          {
            reg_param[npoints] =gamma[p];
            ratio = (gamma[N-M+1]/reg_param[npoints])^(1/(npoints-1));
          }
      }


    if(length(varargin)==2)
      {
        alpharange=varargin
        reg_param[npoints]=alpharange[2];
        ratio=(alpharange[1]/alpharange[2])^(1/(npoints-1));
      }

### % calculate all the regularization parameters
        for(i in seq(from=npoints-1, by=-1, to=1))
          {
            reg_param[i] = ratio*reg_param[i+1];
          }

    if(M>N)
      {
        k=0;
      }
    else
      {
        k=N-M;
      }
    Y=t( solve(X)) ;
    gamma=lambda/MYmu;
    
    ng=length(gamma);
    

### %solve for each solution
    for(i in 1:npoints)
      {
### %evaluate series filter coefficients for this regularization parameter
        f=rep(0, ng);
        for(j in 1:ng)
          {
            f[j]=gamma[j]^2/(gamma[j]^2+reg_param[i]^2);
            if(is.nan(gamma[j]) | is.infinite(gamma[j] ))
              {
                f[j]=1;
              }
            
            if(lambda[j] == 0 &  mu[j] == 0)
              {
                f[j]=0;
              }

### %build the solution

            m[,i]=m[ ,i]+f[j]*(t(U[ ,j+k])%*%d/lambda[j]) %*% Y[,j];
          }

        rho[i]=Mnorm(G%*%m[,i]-d);
        eta[i]=Mnorm(L%*%m[,i]);     
      }

    return(list(rho=rho,eta=eta,reg_param=reg_param,m=m))


  }

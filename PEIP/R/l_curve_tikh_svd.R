l_curve_tikh_svd <-
function(U,s,d,npoints,varargin=NULL)
  {


    eps = .Machine$double.neg.eps
    
###  [rho,eta,reg_param] = l_curve_tikh_svd
### % Initialization. 
    d1 = dim(U);
    m = d1[1]
    n = d1[2]

    p = length(s); 
    
### % compute the projection, and residual error introduced by the projection
    d_proj = t(U) %*% d;
    dr = Mnorm(d)^2 - Mnorm(d_proj)^2;

### %data projections
    d_proj = d_proj[1:p]; 

### %scale series terms by singular values
    d_proj_scale = d_proj/s; 
    
### % initialize storage space
    eta = rep(0, npoints);
    rho = eta;
    reg_param = eta;
    s2 = s^2; 

    if(length(varargin)==0)
      {
### % set the smallest regularization parameter that will be used
        smin_ratio = 16*eps;
        reg_param[npoints] = max(c( s[p],s[1]*smin_ratio) ) ; 

### % ratio so that reg_param(1) will be s(1)
        ratio = (s[1]/reg_param[npoints])^(1/(npoints-1));
      }

    if(length(varargin)==2)
      {
        alpharange=as.vector(varargin);
        reg_param[npoints]=alpharange[2];
        ratio=(alpharange[1]/alpharange[2])^(1/(npoints-1));
      }
    

### % calculate all the regularization parameters
    for(i in seq(from=npoints-1, by=-1, to=1))
      {
        reg_param[i] = ratio %*% reg_param[i+1];
      }

### % determine the fit for each parameter
    for(i in 1:npoints)
      {
### %GSVD filter factors
        f = s2/(s2 + reg_param[i]^2); 
        eta[i] = Mnorm(f*d_proj_scale); 
        rho[i] = Mnorm((1-f)*d_proj); 
      }

### % if we couldn't match the data exactly add the projection induced misfit

    if (m > n && dr > 0)
      {
        rho = sqrt(rho^2 + dr);
      }


    return(list(rho=rho, eta=eta, reg_param=reg_param ) )
  }

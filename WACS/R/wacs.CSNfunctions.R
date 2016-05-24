  ###################################################################
  #
  # These functions are part of WACSgen V1.0 
  # Copyright © 2013,2014,2015, D. Allard, BioSP,
  # and Ronan Trepos MIA-T, INRA
  #
  # This program is free software; you can redistribute it and/or
  # modify it under the terms of the GNU General Public License
  # as published by the Free Software Foundation; either version 2
  # of the License, or (at your option) any later version.
  #
  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details. http://www.gnu.org
  #
  ###################################################################

  ###################################################################
  #
  #  This file contains several functions for computing quantities related
  #  to CSN distributions: density, cumulative probabilities, quantiles, 
  #  marginal densities and marginal parameters. Random draws of CSN
  #  distributions and conditional simulations are also possible
  #' @import mvtnorm  
  ###################################################################
  
  
  
  
  ###################################################################
  #
  #
  # Set of CSN functions for computing the density, cumulative function or quantiles
  # of a univariate RV
  #     CNS(location,scale^2,skew)= CSN_{1,1}(location,scale^2,skew/scale,0,1-skew^2)
  #
  # Note: scale^2 must be thought of as variance, and skew is in [-1,1].
  #
  # ARGUMENT 
  #    z or p: value at which the function is computed
  #    n    : number of random draws
  #
  # VALUE
  #   density, cpf, quantile or random draw
  #
  ###################################################################

  dcsn = function (z, location=0, scale=1, skew=0){  
    # Computes the density of a CNS*(location,scale^2,shape) = CSN_{1,1}(location,scale^2,\delta,1-d^2)
    dz    = 2*dnorm(z,location,scale)*pnorm(skew*(z-location)/scale,0,sqrt(1-skew^2))
  return(dz)
  }

  pcsn = function (z, location=0, scale=1, skew=0){
    # Computes the cumulative probability of a CNS*(location,scale^2,shape) = CSN_{1,1}(location,scale^2,\delta,1-d^2)
    S     = matrix(c(scale^2,-scale^2*skew,-scale^2*skew,1),ncol=2)
    Fz    = 2*pmnorm(c(z,0),c(location,0),S)[1]
  return(Fz)
  }

  qcsn = function(p, location = 0, scale = 1, skew =0){
    # Computes the quantile of a CNS*(location,scale^2,shape) = CSN_{1,1}(location,scale^2,\delta,1-d^2)
    g = function(z){(pcsn(z,location, scale, skew)-p)^2}
    q = uniroot(g,interval=c(qnorm(p/2,location, scale),qnorm(1-p/2,location, scale)))$root
    return(q)
  }

  rcsn = function(n=1, location=0, scale=1, skew=0){
    # random generation of a CNS*(location,scale^2,shape) = CSN_{1,1}(location,scale^2,\delta,1-d^2)
    U.pos = abs(rnorm(n))
    V      = rnorm(n)  
    Z      = location + skew*scale*U.pos + sqrt(1-skew^2)*scale*V
    return(Z)
  }


  ###################################################################
  #
  #
  # Function for computing the density of a  multivariate CNS*(Mu,Sigma,Skew) 
  # or for drawing a random CSN* vector, with
  #     CNS*(Mu,Sigma,Skew)  = CSN_{k,k}(Mu,Sigma,Skew%*%Sigma^{-1/2},0,Id_k-Skew^2) 
  # where
  #     k     dimension of Mu
  #     Mu    vector of location parameter
  #     Sigma matrix of Covariance matrix
  #     Skew  vector of skewness parameters (must be in [-1,1]) 
  #
  # See Flecher et al. 2010
  #
  # ARGUMENT 
  #     Z              : vector at which the density is computed
  #     Mu,Sugma,Skew  : parameters of the CSN* distribution
  #     n              : number of random draws
  #
  #  VALUE
  #   density or random draws
  #
  ###################################################################

  dmcsnstar=function(Z,Mu,Sigma,Skew){
    N = length(Z)
    if( length(Mu)!=N | (dim(Sigma)[1]!=N) | (dim(Sigma)[2]!=N) |  length(Skew)!=N){
      stop ("[dmcsnstar] dimension error")  
    }else{}
    D = diag(Skew)%*%solve(sqrtm(Sigma))
    Delta=diag(N) - diag(Skew)^2
  
    #Computation
    Sigma = round(Sigma,digits=6); Delta = round(Delta,digits=6)
    d = dmnorm(Z,Mu,Sigma)*pmnorm(as.vector(D%*%(Z-Mu)),rep(0,N),Delta)[1]
    d = d*2^(N)
    return(d)
  }

  rmcsnstar=function(n=1,Mu,Sigma,Skew){

    Nv = length(Mu)
    if( (dim(Sigma)[1]!=Nv) | (dim(Sigma)[2]!=Nv) |  length(Skew)!=Nv){ 
      stop ("[rmcsnstar] dimension error")
    }
    z  = matrix(0,n,Nv)
    NonSkew = sqrt(1-Skew^2)
    SS      = sqrtm(Sigma)
    for (i in 1:n){
      u       = abs(rnorm(Nv))
      v       = rnorm(Nv)      
      y       = Skew*u + NonSkew*v
      z[i,]   = Mu + SS%*%y
    }
    return(z)
  }



  ###################################################################
  #
  #
  #  Functions for computing the density of a multivariate CSN(Mu,Sigma,D,Nu,Delta) 
  #  and for drawinf random vectors.
  #  Coded according to Karimi and Mohammadzadeh (2012), Proposition 1(ii)
  #
  # ARGUMENT 
  #     n     number of random draws
  #     y     vector at which the density is computed 
  #     Mu    vector of joint location parameters
  #     Sigma matrix of joint covariance matrix
  #     D     matrix with skewness parameters
  #     Nu    vector of joint location parameters
  #     Delta matrix of joint location parameters
  #     
  #
  # VALUE
  #     value of the density (dmcsn), or n random draws (rmcsn)
  #
  ###################################################################
  rmcsn = function(n,Mu,Sigma,D,Nu,Delta){
  
  # Check dimensions
  if ( (dim(Sigma)[1] != dim(Sigma)[2]) | (dim(Delta)[1] != dim(Delta)[2]) ){ 
    stop ("[rmcsn] covariance matrices must be squared matrices")
  }
  if ( (length(Mu) !=  dim(Sigma)[1]) | (length(Mu) !=  dim(D)[2]) | (dim(Sigma)[1] !=   dim(D)[2]) ){  
    stop ("[rmcsn] dimension error -- vector Mu not compatible with Sigma or D")
  }
  if ( (length(Nu) !=  dim(Delta)[1]) | (length(Nu) != dim(D)[1]) | (dim(D)[1] != dim(Delta)[1] ) ){
    stop ("[rmcsn] dimension error -- vector Nu not compatible with Delta or D")
  }
  
  # Computes parameters of the marginal distribution
  Nv = length(Mu)
  QQ  = Delta + D%*%Sigma%*%t(D)
  FF = Sigma%*%t(D)%*%solve(QQ)
  GG = sqrtm(Sigma-FF%*%D%*%Sigma)
  
  # Does the computation
  OK = FALSE
  add = 0
  niter = 0
  while (!OK){
    niter = niter + 1
    add = add + 0.1
    u.start = Nu + abs(rnorm(length(Nu),0,add))
    u  = rtmvnorm(n,sigma=QQ,lower=Nu,algorithm="gibbs",start.value=u.start,burn.in.samples=200)
    if (!is.na(u[1])) OK = TRUE
    if (niter == 20) {
      u = rmvnorm(n,sigma=QQ)
      u = u + Nu
      OK = TRUE
    }
  }
  v  = matrix(rnorm(length(Mu)*n),nrow=n)
  z  = matrix(rep(Mu,n),nrow=n,byrow=T) + u%*%t(FF) + v%*%GG
  return(z)  
  }
  


  dmcsn=function(W,Mu,Sigma,D,Nu,Delta){
  # Check dimensions
  N = length(W)
  if( (length(Mu)!=N) | (dim(Sigma)[1]!=N) | (dim(Sigma)[2]!=N) | (dim(D)[2]!=N) ){
    stop("[dmcsn] dimension error")
  }
  if( (dim(D)[1]!= dim(Delta)[1] ) | (dim(D)[1]!= dim(Delta)[2] ) ){ 
    stop("[dmcsn] dimension error")
  }
  
  # Does the computation
  d = dmnorm(W,Mu,Sigma)*pmnorm(as.vector(D%*%(W-Mu)),Nu,Delta)[1]
  d = d*2^(-N)
  return(d)
  }



dmcsn.marg=function(y,sel,Mu,Sigma,D,Nu,Delta){
 ####################################################################
 #
 #
 # Computes the marginal density of a CSN distribution with parameters
 #     Mu    vector of joint location parameters
 #     Sigma matrix of joint covariance matrix
 #     D     matrix with skewness parameters
 #     Nu    vector of joint location parameters
 #     Delta matrix of joint location parameters
 #
 #  Coded according to Dominguez-Molina, Gonzales-Farias and Gupta, 2003
 #
 # ARGUMENT 
 #     sel    selection of variables ofr the marginal density
 #       y    vector at which the density is computed 
 # VALUE
 #     value of the density 
 #
 # CALLS
 #     pmnorm and dmnorn from the package mnormt
 ####################################################################
 
  # Check dimensions
  N.y  = length(y)
  N.Mu = length(Mu)
  if (N.y > N.Mu) stop ("[dmcsn.marg] dimension of y must be < dimension of Mu")
  if (N.y != length(sel)) stop ("[dmcsn.marg] dimension of y must be = length(sel)")
  if( (dim(Sigma)[1]!=N.Mu) | (dim(Sigma)[2]!=N.Mu) | (dim(D)[2]!=N.Mu) )    stop("[dmcsn.marg] dimension error between Sigma and D")
  if( (dim(D)[1]!= dim(Delta)[1] ) | (dim(D)[1]!= dim(Delta)[2] ) ) stop("[dmcsn.marg] dimension error between D and Delta")
  

  # Computes parameters of the marginal distribution
  Mu.1       = Mu[sel]
  Sigma.11   = Sigma[sel,sel,drop=FALSE]
  Sigma.12   = Sigma[sel,-sel,drop=FALSE]
  Sigma.22   = Sigma[-sel,-sel,drop=FALSE]
  Sigma.cond = Sigma.22- t(Sigma.12)%*%solve(Sigma.11)%*% Sigma.12
  Sigma.11   = round(Sigma.11,digits=6)
  D.1        = D[,sel,drop=FALSE]
  D.2        = D[,-sel,drop=FALSE]
  #print(" dbg dmcsn.marg")
  #str(Sigma.11)
  D.star     = D.1 + D.2%*%t(Sigma.12)%*%solve(Sigma.11)
  Delta.star = Delta + D.2%*%Sigma.cond%*%t(D.2)
  Delta.star = round(Delta.star,digits=6)
  DD         = Delta.star + D.star%*%Sigma.11%*%t(D.star)  

  # Does the computation
  d  = dmnorm(y,Mu.1,Sigma.11)*pmnorm(as.vector(D.star%*%(y-Mu.1)),Nu,Delta.star)[1]
  d  = d/pmnorm(rep(0,length(Nu)),Nu,DD)[1]
  return(d)
}



extract.csn.marg=function(sel,Mu,Sigma,D,Nu,Delta){
  ###################################################################
  #
  #
  # Extracts the parameters of a marginal density of a CSN distribution with parameters
  #     Mu    vector of joint location parameters
  #     Sigma matrix of joint covariance matrix
  #     D     matrix with skewness parameters
  #     Nu    vector of joint location parameters
  #     Delta matrix of joint location parameters
  #
  #  Coded according to Dominguez-Molina, Gonzales-Farias and Gupta, 2003
  #
  # ARGUMENT 
  #     sel    selection of variables ofr the marginal density
  #       
  # VALUE
  #     A list containing Mu.1, Sigma.11, D.start, Delta.star, DD
  #
  ###################################################################
  
  # Check dimensions
  N.Mu = length(Mu)
  if( (dim(Sigma)[1]!=N.Mu) | (dim(Sigma)[2]!=N.Mu) | (dim(D)[2]!=N.Mu) )    stop("dmcsn.marg: dimension error")
  if( (dim(D)[1]!= dim(Delta)[1] ) | (dim(D)[1]!= dim(Delta)[2] ) ) stop("dmcsn.marg: dimension error")
  
  
  # Computes parameters of the marginal distribution
  Mu.1       = Mu[sel]
  Sigma.11   = Sigma[sel,sel,drop=FALSE]
  Sigma.12   = Sigma[sel,-sel,drop=FALSE]
  Sigma.22   = Sigma[-sel,-sel,drop=FALSE]
  Sigma.cond = Sigma.22- t(Sigma.12)%*%solve(Sigma.11)%*% Sigma.12
  Sigma.11   = round(Sigma.11,digits=6)
  D.1        = D[,sel,drop=FALSE]
  D.2        = D[,-sel,drop=FALSE]
  D.star     = D.1 + D.2%*%t(Sigma.12)%*%solve(Sigma.11)
  Delta.star = Delta + D.2%*%Sigma.cond%*%t(D.2)
  Delta.star = round(Delta.star,digits=6)
  DD         = Delta.star + D.star%*%Sigma.11%*%t(D.star)  
  
  return(list(Mu.1=Mu.1,Sigma.11=Sigma.11,D.star=D.star,Delta.star=Delta.star,DD=DD))
}


rmcsn.cond=function(n,y.cond,sel.cond,Mu,Sigma,D,Nu,Delta){
  ###################################################################
  #
  #
  # performs n random draws of a conditional multivariate CSN distribution:
  #     - length of y (say n.o) must be < length of Mu
  #     - it is assumed that the n.o first coordinates of mu correspond to that of y.
  #     - it is also assumed that the n.o x n.o upper block of Sigma corresponds to that of y
  #     - similar requirement for D  
  # 
  #
  #
  # ARGUMENT 
  #     n        number of random draws
  #     y.cond   vector of conditioning values
  #     sel.cond selection of variables corresponding to y.cond
  #     Mu    vector of location parameters
  #     Sigma matrix of covariance matrix
  #     D     matrix with skewness parameters
  #     Nu    vector of location parameters
  #     Delta matrix of location parameters
  #
  # VALUE
  #     n random draws 
  #
  # CALLS 
  #     rmcsn
  #
  ##########################################################################
  
  # Check dimensions
  N.y  = length(y.cond)
  N.Mu = length(Mu)
  if ( (dim(Sigma)[1] != dim(Sigma)[2]) | (dim(Delta)[1] != dim(Delta)[2]) ) 
    stop ("[rmcsn.cond] covariance matrices must be squared")
  if (any(eigen(Sigma)$values<0)) stop  ("[rmcsn.cond] Sigma be positive definite")
  if (any(eigen(Delta)$values<0)) stop  ("[rmcsn.cond] Sigma be positive definite")
  if (N.y > N.Mu) stop ("[rmcsn.cond] dimension of y must be < dimension of Mu")
  if (N.y != length(sel.cond)) stop ("[rmcsn.cond] dimension of y must be = length(sel)")
  if ( (length(Mu) <= length(y.cond)) | (dim(Sigma)[1] <= length(y.cond)) |  (dim(Delta)[1] <= length(y.cond) ) )
    stop ("[rmcsn.cond] dimension error -- vector y is too long")
  if ( (length(Mu) !=  dim(Sigma)[1]) | (length(Mu) !=  dim(D)[2]) | (dim(Sigma)[1] !=   dim(D)[2]) )   
    stop ("[rmcsn.cond] dimension error -- vector Mu not compatible with Sigma or D")
  if ( (length(Nu) !=  dim(Delta)[1]) | (length(Nu) != dim(D)[1]) | (dim(D)[1] != dim(Delta)[1] ) ) 
    stop ("[rmcsn.cond] dimension error -- vector Nu not compatible with Delata or D")
  
  # Computes parameters of the marginal distribution

  Sigma.11   = Sigma[sel.cond,sel.cond,drop=FALSE]
  Sigma.12   = Sigma[sel.cond,-sel.cond,drop=FALSE]
  Sigma.22   = Sigma[-sel.cond,-sel.cond,drop=FALSE]
  Sigma.cond = Sigma.22 - t(Sigma.12)%*%solve(Sigma.11)%*%Sigma.12
  Sigma.cond = round(Sigma.cond,digits=6)
  
  Mu.cond    = Mu[-sel.cond,drop=FALSE] + t(Sigma.12)%*%solve(Sigma.11)%*%(y.cond-Mu[sel.cond])
  D.1        = D[,sel.cond,drop=FALSE]
  D.2        = D[,-sel.cond,drop=FALSE]
  D.star     = D.1 + D.2%*%t(Sigma.12)%*%solve(Sigma.11)
  Nu.cond    = Nu - D.star%*%(y.cond-Mu[sel.cond])
  
  # Does the simululation
  z = rmcsn(n,as.vector(Mu.cond),Sigma.cond,D.2,as.vector(Nu.cond),Delta)
  return(z)
} 

  
  ###################################################################
  #
  #
  # Functions for computing the 1d and 2d densities of mixtures of CNS
  #     
  #
  # Note: scale^2 must be thought of as variance, and skew is in [-1,1].
  #
  # ARGUMENT 
  #    V, V1, V2:  variables of which the density is computed 
  #    uu       :  value at which the density is computed
  #    Season.par : parameters of the seasons
  #    Prob       : marginal probabilities of the weather states (dry or wet)
  #    vdry       : indicator of dry/wet state 
  # VALUE
  #   density
  #
  ###################################################################  
  
  wacsdensity_1d = function(V,uu,wt,Season.par,Prob,vdry){
    h  = rep(0,length(uu))
    for (w in wt){
      wt.par = Season.par[[paste("W",sep="",w)]]
      Mu     = wt.par$loc
      Skew   = wt.par$skew
      Cov    = wt.par$cov
      Rho    = wt.par$rho
      SS     = Cov
      D      = diag(Skew)%*%sqrtm(solve(SS))
      Delta  = diag(length(Skew)) - diag(Skew^2)
      for (i in 1:length(uu)){
        h[i] = h[i] + Prob[w]*dmcsn.marg(uu[i],V-vdry,Mu,SS,D,rep(0,length(Mu)),Delta)
      }
    }
    h = h/sum(Prob[wt])
    return(h)
  }

  wacsdensity_time = function(V1,uu,wt,Season.par,Prob,vdry){
    g  = matrix(0,length(uu),length(uu))

    for (w in wt){
      wt.par = Season.par[[paste("W",sep="",w)]]
      Nv     = length(wt.par$loc)
      Mu     = c(wt.par$loc,wt.par$loc)
      Skew   = c(wt.par$skew,wt.par$skew)
      Cov    = wt.par$cov
      Rho    = diag(wt.par$rho)
      M      = cbind(Cov,Cov%*%Rho)
      M      = rbind(M,cbind(Rho%*%Cov,Cov))
      Q      = sqrtm(def.pos(solve(M)))
      D      = diag(Skew)%*%Q
      Delta  = diag(length(Skew)) - diag(Skew)%*%diag(Skew)

      for (i in 1:length(uu)){
        for (j in 1:length(uu)){
          g[i,j] = g[i,j] + Prob[w]*dmcsn.marg(c(uu[i],uu[j]),c(V1-vdry,Nv+V1-vdry),Mu,M,D,rep(0,2*Nv),Delta)
        }
      }
    }
    g = g/sum(Prob[wt])
    return(g)
  }
  
  wacsdensity_2d = function(V1,V2,uu,wt,Season.par,Prob,vdry){
    g  = matrix(0,length(uu),length(uu))
    for (w in wt){
      wt.par = Season.par[[paste("W",sep="",w)]]
      Mu     = wt.par$loc
      Skew   = wt.par$skew
      Cov    = wt.par$cov
      Rho    = wt.par$rho
      SS     = Cov
      D      = diag(Skew)%*%sqrtm(solve(SS))
      Delta  = diag(length(Skew)) - diag(Skew^2)
      for (i in 1:length(uu)){
        for (j in 1:length(uu)){
          g[i,j] = g[i,j] + Prob[w]*dmcsn.marg(c(uu[i],uu[j]),c(V1-vdry,V2-vdry),Mu,SS,D,rep(0,length(Mu)),Delta)
        }
      }
    }
    g = g/sum(Prob[wt])
    return(g)
  }
  
  
  
  
  
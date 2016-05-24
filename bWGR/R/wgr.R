wgr = function(y,gen,it=1500,bi=500,th=1,bag=1,
               rp=TRUE,iv=FALSE,pi=0,df=5,R2=0.5,
               eigK=NULL,EigT=0.05,RO=FALSE,verb=FALSE)
{
  
  anyNA = function(x) any(is.na(x))
  if(anyNA(gen)){stop('No missing values allowed in Z')}
  if(anyNA(pi>=1)){stop('Pi must be lower than 1')}
  gen0 = gen
  
  # Polygenic term
  if(!is.null(eigK)){
    pk = sum(eigK$values>EigT)
    U0 = U = eigK$vectors[,1:pk]
    V = eigK$values[1:pk]
    H = h = rep(0,pk)
    Vk = rep(1,pk)
    xxK = rep(bag,pk)
  }  
  
  # Remove missing y's
  if(anyNA(y)){
    mis = which(is.na(y))
    y = y[-mis]
    gen = gen[-mis,]
    if(!is.null(eigK)) U = U[-mis,]
  }
  
  # MCMC settings
  post = seq(bi,it,th)
  mc = length(post)
  
  # Preparing inputs
  X = gen
  n = nrow(gen)
  p = ncol(gen)
  xx = apply(X,2,crossprod)*bag
  b = g = rep(0,p)
  d = rep(1,p)
  mu = mean(y - X%*%b)
  e = y - X%*%g - mu
  Va = 0.0001
  Vb = rep(Va,p)
  Ve = 1
  md = pi
  
  # Priors
  df_prior = df 
  MSx = sum(apply(X,2,var,na.rm=T))
  S_prior = R2*var(y,na.rm=T)*(df_prior+2)/MSx/(1-pi)
  shape_prior  = 1.1
  rate_prior = (shape_prior-1)/S_prior
  S_conj = rgamma(1,p*df_prior/2+shape_prior,sum(1/Vb)/2+rate_prior)
  if(!is.null(eigK)) Sk_prior = R2*var(y,na.rm=T)*(df_prior+2)/pk
  
  # Storing Posterior
  B0 = VA = VE = VP = S = 0
  VB = G = D = B = rep(0,p)
  
  
  #RUN
  if(verb) pb = txtProgressBar(style = 3)
  
  for(i in 1:it){
    
    if(RO){
      ro = sample(1:p)-1
      if(!is.null(eigK)) rok = sample(1:pk)-1
    }
    
    # Resampling
    if(bag!=1) Use = sort(sample(n,n*bag,rp))
    
    # Update polygenic term and regression coefficients
    if(!is.null(eigK)){
      
      Lk = Ve/(Vk*V)
      if(bag!=1){
        
        if(RO){
          update = KMUP2(U[Use,],h,xxK,e[Use],Lk,pk,Ve,0,rok)
        }else{
          update = KMUP(U[Use,],h,xxK,e[Use],Lk,pk,Ve,0)
        }
        
      }else{
        
        if(RO){
          update = KMUP2(U,h,xxK,e,Lk,pk,Ve,0,rok)
        }else{
          update = KMUP(U,h,xxK,e,Lk,pk,Ve,0)
        }      
        
      }
      
      h = update[[1]]
      e = update[[3]]
      
      if(pi>0){
        PI = rbeta(1,10*pi+md+1,10*(1-pi)-md+1)
      }else{
        PI=0
      }
      
      L = Ve/Vb
      if(bag!=1){
        
        if(RO){
          update = KMUP2(X[Use,],b,xx,e,L,p,Ve,PI,ro)
        }else{
          update = KMUP(X[Use,],b,xx,e,L,p,Ve,PI)
        }
        
      }else{
        
        if(RO){
          update = KMUP2(X,b,xx,e,L,p,Ve,PI,ro)
        }else{
          update = KMUP(X,b,xx,e,L,p,Ve,PI)
        }
        
      }
      b = update[[1]]
      d = update[[2]]; if(pi>0) d[is.nan(d)] = 1
      e = update[[3]]
      g = b*d
      md = mean(d)
      
    }else{
      
      # Update regression coefficients without polygene
      if(pi>0){
        PI = rbeta(1,10*pi+md+1,10*(1-pi)-md+1)
      }else{
        PI=0
      }
      
      L = Ve/Vb
      if(bag!=1){
        
        if(RO){
          update = KMUP2(X[Use,],b,xx,e[Use],L,p,Ve,PI,ro)
        }else{
          update = KMUP(X[Use,],b,xx,e[Use],L,p,Ve,PI)
        }
        
      }else{
        
        if(RO){
          update = KMUP2(X,b,xx,e,L,p,Ve,PI,ro) 
        }else{
          update = KMUP(X,b,xx,e,L,p,Ve,PI)
        }
        
      }
      b = update[[1]]
      d = update[[2]]; if(pi>0) d[is.nan(d)] = 1
      e = update[[3]]
      g = b*d
      md = mean(d)
      
    }
    
    # Update genetic variance
    if(iv){    
      Vb = (S_conj+b^2)/rchisq(p,df_prior + 1)
      S_conj = rgamma(1,p*df_prior/2+shape_prior,sum(1/Vb)/2+rate_prior)
    }else{
      Va = (sum(b^2)+S_prior)/rchisq(1,df_prior + p)
      Vb = rep(Va,p)
    }
    if(!is.null(eigK)){
      Vp = (sum(h^2/V)+Sk_prior)/rchisq(1,df_prior + pk)
      Vk = rep(Vp,pk)
    }
    
    
    # Residual variance
    Ve = (sum(e*e)+S_prior)/rchisq(1,n*bag + df_prior)
    
    # Intercept
    if(!is.null(eigK)){e = y-X%*%g-U%*%h}else{e = y-X%*%g}
    mu = rnorm(1,mean(e),(Ve+1E-5)/n)
    if(is.na(mu)||is.nan(mu)) mu = mean(y,na.rm=T)
    e = e-mu
    
    # Save posterior
    if(i%in%post){
      B0 = B0+mu
      B = B+b
      D = D+d
      G = G+g
      VE = VE+Ve
      if(iv){VB = VB+Vb}else{VA = VA+Va}
      if(!is.null(eigK)){H = H+h; VP = VP+Vp}
    }
    
    if(verb) setTxtProgressBar(pb, i/it)
  }
  
  if(verb) close(pb)
  
  # Posterior mean
  B0 = B0/mc
  B = B/mc
  D = D/mc
  G = G/mc
  VE = VE/mc
  if(iv){VB = VB/mc}else{VA = VA/mc}
  if(!is.null(eigK)){H = H/mc; VP = VP/mc}
  
  # Fitted values
  if(!is.null(eigK)){
    poly = U0 %*% H
    HAT = B0 + gen0 %*% G + poly
  }else{
    HAT = B0 + gen0 %*% G
  }
  
  
  # Output
  if(!is.null(eigK)){
    if(iv){
      final = list('mu'=B0,'b'=B,'Vb'=VB,'g'=G,'d'=D,'Ve'=VE,'hat'=HAT,'u'=poly,'Vk'=VP)
    }else{
      final = list('mu'=B0,'b'=B,'Vb'=VA,'g'=G,'d'=D,'Ve'=VE,'hat'=HAT,'u'=poly,'Vk'=VP)
    }
  }else{
    if(iv){
      final = list('mu'=B0,'b'=B,'Vb'=VB,'g'=G,'d'=D,'Ve'=VE,'hat'=HAT)
    }else{
      final = list('mu'=B0,'b'=B,'Vb'=VA,'g'=G,'d'=D,'Ve'=VE,'hat'=HAT)
    }
  }
  
  
  
  return(final)
}
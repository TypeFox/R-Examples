lk_obs_manifest <- function(par,Y,Xd,indn,lev,k,sup,G2,IPI,mod,outp=FALSE){

#        [lk,U,s] = lk_obs(par,Y,X,Xd,indn,lev,k,sup,G2,IPI,mod=0)
#
# INPUT:
# Y:    matrix of response variables
# D:    design matrix for the covariates
# X:    matrix of covariates
# lev:  vector containing the number of levels of each variable
# type: vector with elements 'l', 'g', 'c', 'r' indicating the type of logit
# k:    number of latent classes
# par:  vector of regression parameters
# La:   starting value for initial probability vector
# PI:   starting value for transition matrix
#
# OUTPUT:
# lk:   log-likelihood
# s:    score of the observed log-likelihood

  # preliminaries
  n = dim(Y)[2]; TT = dim(Y)[3]; nd = dim(Xd)[3]
  np = lev-1
  q = length(sup)
  Gm = cbind(-rep(1,lev-1,1),diag(lev-1))
  Hm = rbind(rep(0,lev-1),diag(lev-1))
  GHt = t(Gm)%*%t(Hm)
  lm = c(1,rep(0,lev-1))
  Lm = rbind(rep(0,lev-1),diag(lev-1))-rbind(diag(lev-1),rep(0,lev-1))
  INDN = vector("list",nd)
  for(jd in 1:nd) INDN[[jd]]$ind = which(as.vector(indn)==jd)
  # invert parametrization
  out = trans_par(par,lev,k,sup,G2,IPI,mod)
  la = out$la; PI = out$PI; rho = out$rho; si = out$si; par = out$par; lrho = out$lrho; tau = out$tau 
  # AR probabilities
  wei = dnorm(sup); wei = wei/sum(wei);
  las = as.vector(la%x%wei)
  if(q>1){
    SUP = sup%o%rep(1,q)
    WEI = matrix(0,k*q,k*q)
    for(j in 1:k){
      ind = (j-1)*q+(1:q)
      Wei = dnorm(t(SUP),rho[j]*SUP,sqrt(1-rho[j]^2))
      Wei = Wei/rowSums(Wei)
  	  WEI[,ind] = matrix(1,k,1)%x%Wei
    }
  	PIs = (PI%x%matrix(1,q,q))*WEI
  }else PIs = PI
  LLm1 = array(0,c(dim(Lm)[2],dim(Lm)[1],nd))
  for (re in 1:nd) LLm1[,,re] = t(Lm)
  # compute initial log-likelihood
  I = diag(np);
  one = matrix(1,np,1);
  Pio = array(0,c(n,k*q,TT))
  par0 = par[1:(lev-1+k)]
  Eta01 = prod_array(Xd,par[(lev+k):length(par)]); j = 0
  for(c in 1:k){
    u = matrix(0,1,k); u[c] = 1; u = u[-1]
    D0 = cbind(I,matrix(u,nrow=1)%x%one)
    for(d in 1:q){
      j = j+1;
      D = cbind(D0,sup[d]*one); agg = D%*%par0
      Eta1 = Eta01+agg%*%matrix(1,1,nd)
      Qv1 = expit(Eta1); Qv1 = pmin(pmax(Qv1,10^-100),1-10^-100)
      Pv1 = lm%*%matrix(1,1,nd)+Lm%*%Qv1; Pv1 = pmin(pmax(Pv1,10^-100),1-10^-100);
      for(t in 1:TT) Pio[,j,t] = colSums(Y[,,t]*Pv1[,indn[,t]])
    }
  }
  # apply recursion
  if(k==1) PIs = as.matrix(PIs)
  Q = rec1(Pio,las,PIs)
  if(q*k==1) pim = Q[,,TT] else pim = rowSums(Q[,,TT])
  lk = sum(log(pim))
  # compute U
  if(outp){
  	out = rec3(Q,PIs,Pio,pim)
  	U = out$U; V = out$V
  }else{U=NULL; V=NULL}
  # compute derivative
  s1 = NULL;
  
  if(outp){
  	
    if(k>1){
      if(mod==0){
      	out = stationary(tau,k,G2,IPI)
        d21 = out$d0; d22 = out$d1 
        Mar = diag(k)%x%matrix(1,1,q)
        u1 = Mar%*%rowSums(U[,,1])
        V1 = Mar%*%V%*%t(Mar)
        s2 = d22%*%(u1/la)+d21%*%(V1[IPI]/PI[IPI])
      }
      if(mod==1){
        Mar = diag(k)%x%matrix(1,1,q)
        u1 = Mar%*%rowSums(U[,,1])
        la = as.vector(la)
        Om = diag(la)-la%o%la
        d22 = t(G2)%*%Om
        s2 = d22%*%(u1/la)
      }
    }else{
      s2 = NULL
    }
    # with respect to logit(rho)
    if(q>1){
      s3 = matrix(0,k,1)
      DD1 = matrix(0,k*q,k*q)
      for(j in 1:k){
        ind = (j-1)%*%q+(1:q)
        D1 = (rho[j]+SUP*(t(SUP)-rho[j]*SUP))/(1-rho[j]^2)-(t(SUP)-rho[j]*SUP)^2*rho[j]/(1-rho[j]^2)^2;
        Wei = dnorm(t(SUP),rho[j]*SUP,sqrt(1-rho[j]^2))
        D1 = D1-(rowSums(Wei*D1)/rowSums(Wei))%*%matrix(1,1,q)
        DD1[,ind] = matrix(1,k,1)%x%D1
        for(j1 in 1:k){
          ind1 = (j1-1)%*%q+(1:q);
          s3[j] = s3[j]+2*exp(lrho[j])/(1+exp(lrho[j]))^2*sum(sum(V[ind1,ind]*D1));
        }
      }
    }else{
      s3 = NULL
    }
    # E-step
    U = aperm(U,c(2,1,3))
    s4 = 0; j = 0
    for(c in 1:k){
      u = matrix(0,1,k); u[c] = 1; u = u[-1]
      D0 = cbind(I,matrix(u,nrow=1)%x%one)
      for(d in 1:q){
        j = j+1;
        D = cbind(D0,sup[d]*one); agg = D%*%par0
        Eta1 = Eta01+agg%*%matrix(1,1,nd)
        Qv1 = expit(Eta1); Qv1 = pmin(pmax(Qv1,10^-100),1-10^-100);
        Pit1 = lm%*%matrix(1,1,nd)+Lm%*%Qv1; Pit1 = pmin(pmax(Pit1,10^-100),1-10^-100);
        QQv1 = Qv1*(1-Qv1);
        DPv1 = 1/Pit1;
        RRtc1 = array(0,c(np,lev,nd))
        for(j1 in 1:np) for(j2 in 1:lev) RRtc1[j1,j2,] = QQv1[j1,]*DPv1[j2,]
        RRtc1 = RRtc1*LLm1
        XXRi1 = array(0,c(dim(D)[1],dim(D)[2]+dim(Xd)[2],nd))
        for(h2 in 1:nd){
        		if(lev==2) XXRi1[,,h2] = c(D, Xd[,, h2])
        		else XXRi1[,,h2] = cbind(D,Xd[,,h2])        	
        } 
        XXRi1 = aperm(XXRi1,c(2,1,3)) 
        pc = U[,j,]; pc = as.vector(pc)
        nt = dim(Y)[1]
        YGP = matrix(Y,nt,n*TT)-Pit1[,as.vector(indn)]
        
        for(jd in 1:nd){
          ind = INDN[[jd]]$ind
          pci = pc[ind]
          if(lev==2) XRi = (XXRi1[, , jd] %o% RRtc1[, , jd]) %*% GHt
          else XRi = (XXRi1[,,jd]%*%RRtc1[,,jd])%*%GHt
          if(length(ind)==1){
            s4 = s4+XRi%*%(YGP[,ind]*pci)
          }else{
            s4 = s4+XRi%*%(YGP[,ind]%*%pci)
          }
        }       
      }
    }
   
    s = c(s1,s2,s3,s4)
  }else{
  	s=NULL
  }  
  out = list(lk=lk,U=U,s=s)
  out
}


est_lm_cov_manifest <- function(S,X,k,q,mod,tol=10^-8,maxit=1000,start=0,mu=NULL,al=NULL,
                                  be=NULL,si=NULL,rho=NULL,la=NULL,PI=NULL,output=FALSE,out_se=FALSE){

#
# Fit the model of Bacci, Bartolucci and Pennoni (2014) with global logits
# (when mod = 1)
#
# INPUT:
# S:     array of available configurations (n x TT x r)
# X:     array (n x TT x nc) of covariates with eventually includes lagged
#        response
# k:     number of latent states
# q:     number of support points for AR
# mod:   model (0 = LM with stationary transition, 1 = finite mixture)
# tol:   tolerance for the convergence (optional) and tolerance of conditional probability
#        if tol>1 then return
# start: equal to 1 for random starting values (optional)
# mu:    starting value for mu (optional)
# al:    starting value for al (optional)
# be:    starting value for be (optional)
# si:    starting value for si (optional)
# rho:   starting value for rho (optional)
# la:    starting value for la (optional)
# PI:    starting value for PI (optional)
# output:  TRUE to return additional output
#out_se:   TRUE for computing information
#
# OUTPUT:
#        mu:    vector of cuptpoints
#        al:    support points for the latent states
#        be:    estimate of the vector of regression parameters
#        si:    sigma of the AR process
#        rho:   parameter vector for AR
#        la:    vector of initial probabilities
#        PI:    transition matrix
# lk:    maximum log-likelihood
# np:   number of parameters
# aic:   AIC index
# bic:   BIC index
# sebe:  standard errors for the regression parameters be
# selrho: standard errors for logit type transformation of rho
# PRED0: prediction of latent state
# PRED1: prediction of the overall latent effect

# preliminaries


# *** organize response matrix ***
lev = max(S)+1
if(min(S)>0){
  		cat("|------------------- WARNING -------------------|\n")
  		cat("|The first response category must be coded as 0 |\n")
  		cat("|-----------------------------------------------|\n")	
 } 
nt = prod(lev)
n = nrow(S); TT = ncol(S)
if(is.array(S)){
	r = dim(S)[3]
	if(!is.na(r) & r>1) warning("multivariate data are not allowed; only the first response variable is considered")
	S = matrix(S,n,TT)
} 
if(is.data.frame(S)) warning("Data frame not allowed for S")
if(n!= dim(X)[1]) stop("dimension mismatch between S and X")

Y0 = S+1
S = array(0,c(nt,n,TT))
for(i in 1:n) for(t in 1:TT){
   ind = Y0[i,t]
   S[ind,i,t] = 1
}
if(is.matrix(X)) X = array(X,c(n,TT,1))
nc = dim(X)[3]
ne = lev-1

XX = X
X = array(0,c(ne,nc,n,TT))
for(i in 1:n) for(t in 1:TT){
   if(lev==2) X[,,i,t] = XX[i, t, ]
   else X[,,i,t] = rep(1,ne)%o%XX[i,t,]
}

opt = list(TolFun=10^-6,TolX=10^-6)
opt1 = list(TolFun=10^-6,Display="iter")
out = marg_param(lev,"g")
Cm = out$C; Mm = out$M
Gm = cbind(-rep(1,lev-1),diag(lev-1))
Hm = rbind(rep(0,lev-1),diag(lev-1))
GHt = t(Gm)%*%t(Hm)
lm = c(1,rep(0,lev-1))
Lm = rbind(rep(0,lev-1),diag(lev-1))-rbind(diag(lev-1),rep(0,lev-1))

if(q==1) sup = 0 else{
 # lim = min(sqrt(q),5);
   lim = 5;
  sup = seq(-lim,lim,2*lim/(q-1))
}
Mar = diag(k)%x%matrix(1,1,q) 
if(mod==0 & q>1){
	cat("|------------------------ WARNING --------------------|\n")
    cat("| When mod=0, q has to be set equal to 1              |\n")   
    cat("|-----------------------------------------------------|\n\n")
}  
G2 = NULL; H2 = NULL; IPI = NULL
if(k>1){
  if(mod==0){ 
    for(c in 1:k){
      G2c = diag(k)[,-c]   
      H2c = diag(k)[-c,]
      if (k==2) H2c[c]=-1 else H2c[,c]= -1
      if(is.null(G2)) G2 = G2c else if(k==2) G2 = blkdiag(matrix(G2,ncol=1),matrix(G2c,ncol=1)) else G2 = blkdiag(G2,G2c) 
      if(is.null(H2)) H2 = H2c else if(k==2) H2 = blkdiag(matrix(H2,nrow=1),matrix(H2c,nrow=1))else H2 = blkdiag(H2,H2c)
      IPI = c(IPI,c+seq(0,k*(k-1),k))
    }
  }
  if(mod==1){ 
    G2 = diag(k)[,-1]; if(k==2) G2 = matrix(G2,ncol=1)
    H2 = diag(k); H2[,1] = -1; H2 = H2[-1,]
  }
}
# starting values
mu_inp=mu
if(is.null(mu)){
  Pim = apply(S,c(1,2),sum)+0.05*TT; Eta = Cm%*%log(Mm%*%Pim)
  Eta = Eta%x%matrix(1,1,TT)
  eta = as.vector(Eta)
  Z = matrix(aperm(X,c(1,4,3,2)),n*ne*TT,dim(X)[2])
  Z = cbind(matrix(1,n*TT,1)%x%diag(ne),Z)
  par = ginv(t(Z)%*%Z)%*%t(Z)%*%eta
  mu = par[1:ne]; par = par[-(1:ne)]; be = par
  if(k==1) al = NULL else{
    if(start==1) al = rnorm(k)*k else al = seq(-k,k,2*k/(k-1))
    mu = mu+al[1]
    al = al[-1]-al[1]
  }
  if(k==1) PI = 1 else{
  	if(mod==0){
  		PI = matrix(1,k,k)+9*diag(k) 
  		PI = diag(1/rowSums(PI))%*%PI
  	}
    if(mod==1) PI = diag(k)
  }
  if(start==1){
    la = matrix(runif(k),k,1); la = la/sum(la)
    rho = 2*matrix(runif(k),k,1)-1
    si = runif(1)*5
  }else{
    la = matrix(1,k,1)/k
    rho = matrix(0,k,1)
    si = 3
  }
}  
if(start==2){
	if(is.null(mu_inp)) stop("initial value of the cut-points (mu) must be given in input")
	mu=mu_inp
	if(is.null(be)) stop("initial value of the regression parameters (be) must be given in input")
	be=be
	if(is.null(al)) stop("initial value of the support points (al) must be given in input")
	
	mu = mu+al[1]
	al=al[-1]-al[1]
	if(is.null(la)) stop("initial value of the initial probabilities (la) must be given in input")
	la=la
	if(is.null(PI)) stop("initial value of the transition probabilities (PI) must be given in input")
	PI=PI
	if(mod==0){
		rho = matrix(0,k,1)
		si = 0
	}
	if(mod==1){
		rho = rho
		si = si
	}
	
}

par = c(mu,al,si,be)

if(k==1) tau = NULL else{
	if(mod==0) tau = H2%*%log(PI[IPI]) else tau = H2%*%log(la)
}
wei = dnorm(sup); wei = wei/sum(wei)
las = as.vector(la%x%wei)
lrho = (rho+1)/2
lrho = log(lrho/(1-lrho))
SUP = sup%o%rep(1,q)
WEI = matrix(0,k*q,k*q)
for(j in 1:k){
  ind = (j-1)*q+(1:q);
  Wei = dnorm(t(SUP),rho[j]*SUP,sqrt(1-rho[j]^2))
  Wei = Wei/rowSums(Wei)
  WEI[,ind] = matrix(1,k,1)%x%Wei
}
PIs = (PI%x%matrix(1,q,q))*WEI

t0 = proc.time()

# to do in Fortran
# find non-redundant X configurations (may be very slow)
# Xd = array(X,c(ne,nc,n*TT))
# indn = matrix(1:(n*TT),n,TT)
# for(jd in 1:nd) INDN[[jd]]$ind = jd
X1 = matrix(X,ne*nc,n*TT)
out1 = t(unique(t(X1)))
nd = ncol(out1)
indn = rep(0,n*TT)
INDN = vector("list",nd)
tmp = ne*nc
for(jd in 1:nd){
	ind = which(colSums(X1 == out1[,jd])==tmp)
	indn[ind] = jd
	INDN[[jd]]$ind = ind
}
indn = matrix(indn,n,TT)
Xd = array(out1,c(ne,nc,nd))
#for(jd in 1:nd) INDN[[jd]]$ind = which(indn==jd)
cat(c("n. distinct covariate conf. = ",nd))
#Xd = Xd[,,1:nd]  #to allow for one covariate
LLm1 = array(t(Lm),c(ncol(Lm),nrow(Lm),nd))
# alternate between EM and NR
itg = 0; cont = 1;

while(cont && itg<5){
	
  cont = 0; itg = itg+1;
  # compute initial log-likelihood
  I = diag(ne)
  one = matrix(1,ne,1)
  Pio = array(0,c(n,k*q,TT))
  par0 = par[1:(lev-1+k)]
  Eta01 = prod_array(Xd,par[(lev+k):length(par)]); j = 0
  for(c in 1:k){
    u = matrix(0,1,k); u[c] = 1; u = u[-1]
    D0 = cbind(I,matrix(u,nrow=1)%x%one)
    for(d in 1:q){
      j = j+1;
      D = cbind(D0,sup[d]*one); agg = D%*%par0 
      Eta1 = Eta01+agg%*%rep(1,nd)  
      Qv1 = expit(Eta1); Qv1 = pmin(pmax(Qv1,10^-100),1-10^-100)
      Pv1 = lm%o%rep(1,nd)+Lm%*%Qv1; Pv1 = pmin(pmax(Pv1,10^-100),1-10^-100)
      for(t in 1:TT) Pio[,j,t] = colSums(S[,,t]*Pv1[,indn[,t]])
    }
  }
  Q = rec1(Pio,las,PIs)
  if(q*k==1) pim = Q[,,TT] else pim = rowSums(Q[,,TT])
  lk = sum(log(pim))
  if(tol>1){
  	est = NULL; return
  }
  # start EM
  t0 = proc.time()[1]; tdisp = 5;
  cat("\nEM step:\n")
  if(mod==0){
  	cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat("      k     |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat(sprintf("%11g",c(k,start,0,lk)),"\n",sep=" | ")
  }else{
    cat("----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n");
    cat("     k    |  n. sup   |  max(rho) |  sigma    |  step     |   lk      |  lk-lko   |discrepancy|\n");
    cat("----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n");
  	cat(sprintf("%9g",c(k,q,max(rho),si,0,lk)),"\n",sep=" | ")
  	}
  	
  # iterate until convergence
  it = 0; lko = lk-10^10; dis = 0
  #	while((dis>tol || lk-lko>tol || it ==0) && it <5000){
    while((lk-lko)/abs(lk)>tol & it<maxit){
  	it = it+1
    lko = lk; paro = par; tauo = tau; lrhoo = lrho
    # E-step
    out = rec3(Q,PIs,Pio,pim)
    U = out$U; V = out$V
    # M-step: latent parameters
    if(k>1){
      if(mod==0){
        u1 = Mar%*%rowSums(U[,,1])
        V1 = Mar%*%V%*%t(Mar)
        out = optim(tau,lk_sta,gr=NULL,as.vector(u1),V1,G2,outl=FALSE,method = "BFGS") 
      	tau = out$par
      	out = lk_sta(tau,as.vector(u1),V1,G2,outl=TRUE)			
        flk = out$flk; la = out$la; PI = out$PI
      }
      if(mod==1){
        u1 = Mar%*%rowSums(U[,,1])
        la = u1/n;
        tau = H2%*%log(la)
      }
    }
    if(q>1){
      for(c in 1:k){
        Vc = matrix(0,q,q);
        ind = (c-1)*q+(1:q);
        for(d in 1:k){
          ind1 = (d-1)*q+(1:q)
          Vc = Vc+V[ind1,ind]
        }
        out = optim(lrho[c],lk_ar_rho,gr=NULL,SUP,Vc,outp=FALSE,method = "BFGS")
        lrho[c] = out$par
        out = lk_ar_rho(lrho[c],SUP,Vc,outp=TRUE)
        flk = out$flk; Wei = out$Wei; rho[c] = out$rho 
        WEI[,ind] = matrix(1,k,1)%x%Wei
      }
    }
    las = as.vector(la%x%wei); PIs = (PI%x%matrix(1,q,q))*WEI
    # M-step: regression parameters
    U = aperm(U,c(2,1,3))
    s = 0; FF = 0; j = 0
    for(c in 1:k){
      u = matrix(0,1,k); u[c] = 1; u = u[-1]
      D0 = cbind(I,t(as.matrix(u))%x%one)
      for(d in 1:q){
        j = j+1
        D = cbind(D0,sup[d]*one); agg = as.vector(D%*%par0)
        Eta1 = Eta01+agg%o%rep(1,nd)
        Qv1 = expit(Eta1); Qv1 = pmin(pmax(Qv1,10^-100),1-10^-100)
        Pit1 = lm%o%rep(1,nd)+Lm%*%Qv1; Pit1 = pmin(pmax(Pit1,10^-100),1-10^-100)
        QQv1 = Qv1*(1-Qv1)
        DPv1 = 1/Pit1
        RRtc1 = array(0,c(ne,lev,nd))
        for(j1 in 1:ne) for(j2 in 1:lev) RRtc1[j1,j2,] = QQv1[j1,]*DPv1[j2,]
        RRtc1 = RRtc1*LLm1
        XXRi1 = array(0,c(dim(D)[1],dim(D)[2]+dim(Xd)[2],nd))
        for(h2 in 1:nd){
        	if(lev==2) XXRi1[,,h2] = c(D, Xd[,, h2])
        	else XXRi1[,,h2] = cbind(D,Xd[,,h2])
        }
        XXRi1 = aperm(XXRi1,c(2,1,3))
        pc = U[,j,]; pc = as.vector(pc)
        nt = dim(S)[1]
        YGP = matrix(S,nt,n*TT)-Pit1[,as.vector(indn)]
        Om = array(0,c(lev,lev,nd))
        for(r1 in 1:lev) for(r2 in 1:lev){
            if(r2==r1){
              Om[r1,r2,] = Pit1[r1,]-Pit1[r1,]*Pit1[r2,]
            }else{
              Om[r1,r2,] = -Pit1[r1,]*Pit1[r2,]
            }
        }
        for(jd in 1:nd){
          ind = INDN[[jd]]$ind
          pci = pc[ind]
          if(lev==2) XRi = (XXRi1[, , jd] %o% RRtc1[, , jd]) %*% GHt
          else XRi = (XXRi1[,,jd]%*%RRtc1[,,jd])%*%GHt
          if(length(ind)==1){
            s = s+XRi%*%(YGP[,ind]*pci)
          }else{
            s = s+XRi%*%(YGP[,ind]%*%pci)
          }
          FF = FF+sum(pci)*(XRi%*%Om[,,jd])%*%t(XRi)
        }
      }
    }
    # update parameter vector
    dpar = ginv(FF)%*%s
    mdpar = max(abs(dpar))
    if(mdpar>1) dpar = dpar/mdpar
    par = par+dpar
    si = par[ne+k]
    # compute new log-likelihood
    par0 = par[1:(lev-1+k)]
    Eta01 = prod_array(Xd,par[(lev+k):length(par)]); j = 0
    for(c in 1:k){
      u = matrix(0,1,k); u[c] = 1; u = u[-1]
      D0 = cbind(I,t(as.matrix(u))%x%one)
      for(d in 1:q){
        j = j+1;
        D = cbind(D0,sup[d]*one); agg = as.vector(D%*%par0)
        Eta1 = Eta01+agg%o%rep(1,nd)
        Qv1 = expit(Eta1); Qv1 = pmin(pmax(Qv1,10^-100),1-10^-100);
        Pv1 = lm%o%rep(1,nd)+Lm%*%Qv1; Pv1 = pmin(pmax(Pv1,10^-100),1-10^-100);
        for(t in 1:TT) Pio[,j,t] = colSums(S[,,t]*Pv1[,indn[,t]])
      }
    }
    Q = rec1(Pio,las,PIs)
    if(k*q==1) pim = Q[,,TT] else pim = rowSums(Q[,,TT])
    lk = sum(log(pim))
    # display results
    dis = max(abs(c(par-paro,tau-tauo,lrho-lrhoo)))
#    if((proc.time()[1]-t0)>tdisp){
	#Display output
	if(it/10 == floor(it/10)){
		if(mod==0){
			cat(sprintf("%11g",c(k,start,it,lk,lk-lko,dis)),"\n",sep=" | ")
		}else{
    	cat(sprintf("%9g",c(k,q,max(rho),si,it,lk,lk-lko,round(dis,7))),"\n",sep=" | ")
    	}
    }
#      tdisp = etime(clock,t0)+5
#    }
  }
  # Newton-Rapshon
  par1 = NULL;
  if(k>1) par1 = tau
  if(q>1) par1 = c(par1,lrho)
  par1 = c(par1,par)
  # NR algorithm
  if(mod==1){
 	cat("--------------------------------------------------------------------------------------\n");
 	cat("\nNR step:\n")
    cat("----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n");
    cat("     k    |  n. sup   |  max(rho) |  sigma    |  step     |   lk      |  lk-lko   |discrepancy|\n");
    cat("----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n");
   }
  it = 0; t0 = proc.time()[1]
  while(abs(lk-lko)>10^-5 && it<100 && mod==1){
    lko = lk
    it = it+1
    out = lk_obs_manifest(par1,S,Xd,indn,lev,k,sup,G2,IPI,mod,outp=TRUE)
    nx = length(par1)
    d0 = out$s
    ny = length(d0)
    D = matrix(0,nx,ny)
    for (i in 1:nx){
    	o = matrix(0,nx,1); o[i] = 10^-6
    	out = lk_obs_manifest(par1+o,S,Xd,indn,lev,k,sup,G2,IPI,mod,outp=TRUE)
    	d1 = out$s
  		d = (d1-d0)/10^-6
  		D[i,] = t(d)
    }
    s1 = d0
    J1 = D
    J1 = -(J1+t(J1))/2
 
 	print(c("rcond of information = ",toString(round(rcond(J1),3))))
    dpar1 = ginv(J1)%*%s1;
    mdpar1 = max(abs(dpar1));
    if(mdpar1>0.5) dpar1 = dpar1/mdpar1*0.5
    par1o = par1;
    par1 = par1+dpar1;
    lktmp = lk_obs_manifest(par1,S,Xd,indn,lev,k,sup,G2,IPI,mod)
    lk = lktmp$lk
    cont = 0;
    while(lk<(lko-10^-6)){
      cont = 1;
      dpar1 = dpar1/2;
      par1 = par1o+dpar1;
      lktmp = lk_obs_manifest(par1,S,Xd,indn,lev,k,sup,G2,IPI,mod)
      lk = lktmp$lk
      print(c('halved step',toString(lk-lko)))
    }
    out = trans_par(par1,lev,k,sup,G2,IPI,mod)
    la = out$la; PI = out$PI; rho = out$rho; si = out$si; par = out$par; 
    lrho = out$lrho; tau = out$tau 
    cat(sprintf("%9g",c(k,q,max(rho),si,it,lk,round(lk-lko,7),round(max(abs(par1-par1o)),7))),"\n",sep=" | ")
    #print(c(k,q,max(rho),si,it,lk,lk-lko,max(abs(par1-par1o)),(proc.time()[1]-t0)/it))
   }
  
  if(cont){    
  	las = as.vector(la%x%wei); PIs = (PI%x%matrix(1,q,q))*WEI 
    #las = as.vector(kron(la,wei))
    WEI = matrix(0,k*q,k*q);
    
    for(j in 1:k){
      ind = (j-1)*q+(1:q);
      Wei = dnorm(t(SUP),rho[j]*SUP,sqrt(1-rho[j]^2))
      Wei = Wei/rowSums(Wei)
  	  WEI[,ind] = matrix(1,k,1)%x%Wei
  }
  
    PIs = (PI%x%matrix(1,q,q))*WEI
    tol = tol/2;
  }   
}
# separate parameters and compute aic and bic
mu = par[1:ne]
al = 0
if(k>1) al = c(al,par[(ne+1):(ne+k-1)])
mu = mu+al%*%la
al = al-al%*%la
be = par[(ne+k+1):length(par)]
if(mod==0) np = k*(k-1)
if(mod==1) np = k-1
np = np + (ne+(k-1)+nc) + ((k+1)*(q>1))
if(q==1){
	si=NULL; rho = NULL
}

# compute aic, bic and prediction of latent structure
aic = -2*lk+2*(np)
bic = -2*lk+log(n)*(np)


# prediction
if(output){
  out = lk_obs_manifest(par1,S,Xd,indn,lev,k,sup,G2,IPI,mod,outp=TRUE)
  lk = out$lk; U = out$U
  sup1 = t(Mar)%*%al
 # if(q>1) sup1 = sup1+kron(ones(k,1),sup*si)
  if(q>1) sup1 = sup1+matrix(1,k,1)%x%(sup*si)
  PRED0 = array(0,c(n,k,TT)); PRED1 = matrix(0,n,TT)
  for(t in 1:TT){
    PRED0[,,t] = U[,,t]%*%t(Mar)
    PRED1[,t] = U[,,t]%*%sup1
  }
}
# standard errors
if(out_se){ 					
  out = lk_obs_manifest(par1,S,Xd,indn,lev,k,sup,G2,IPI,mod,outp=TRUE)
  nx = length(par1)
  d0 = out$s
  ny = length(d0)
  D = matrix(0,nx,ny)
  for (i in 1:nx){
    	o = matrix(0,nx,1); o[i] = 10^-6
    	out = lk_obs_manifest(par1+o,S,Xd,indn,lev,k,sup,G2,IPI,mod,outp=TRUE)
    	d1 = out$s
  		d = (d1-d0)/10^-6
  		D[i,] = t(d)
  }
  s1 = d0
  J1 = D
  J1 = -(J1+t(J1))/2 
  print(c("rcond of information = ",rcond(J1)))
  se1 = sqrt(diag(ginv(J1)))
  if(k>1){
  	if(mod==0) se1 = se1[-(1:(k*(k-1)))]
  	if(mod==1) se1 = se1[-(1:(k-1))]
  }
  if(q==1) lrho = NULL
  if(q>1){
  	lrho = se1[1:k]
  	se1 = se1[-(1:k)]
  }
  #se = list(lrho=lrho, be = se1[(ne+k+1):length(se1)])
  selrho = lrho
  sebe = se1[(ne+k+1):length(se1)]
}

out = list(mu=mu,al=al,be=be,si=si,rho=rho,la=la,PI=PI,lk=lk,np=np,aic=aic,bic=bic,call=match.call())
if(out_se){
	out$selrho=selrho
	out$sebe = sebe
	out$J1=J1
}
if(output){
	out$Phi = Pio
	out$PRED0 = PRED0
	out$PRED1 = PRED1
}
 if(mod==0){	
 	cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
 }else{
 	cat("----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n");
 	}	
 	 class(out)="LMmanifest"
	return(out)
}
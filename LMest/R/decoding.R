decoding <- function(est,Y,X1=NULL,X2=NULL,fort=TRUE){
	
	# Provide local decoding on the basis of the output of
	# est_lm_basic, est_lm_cov_latent, est_lm_cov_manifest (in this case the covariates are only in X1)
	# and est_lm_mixed
	
	# est = output from one of these function
	# Y = vector of responses or matrix of responses for which having local deconding
	#
	# Ul = matrix of local deconding corresponding to each row of Y
	# Ug = matrix of global deconding corresponding to each row of Y
	
# est_lm_basic	   
	if(class(est)=="LMbasic"){
		miss = any(is.na(Y))
		if(miss){
        	R = 1 * (!is.na(Y))
         	Y[is.na(Y)] = 0
		}else{
			R = NULL
		}
		if(dim(est$Psi)[3]==1){
			if(is.vector(Y)){
				Y = t(Y)
				if(miss) R = t(R)
			} 		
			n = nrow(Y); TT = ncol(Y)
		}else{
			if(is.matrix(Y)){ 
				Y = array(Y,c(1,dim(Y)))
				if(miss) R = array(R,c(1,dim(R)))
			}
			n = dim(Y)[1]; TT = dim(Y)[2]; r = dim(Y)[3]
		}
		piv = est$piv; Pi = est$Pi; Psi = est$Psi
		k = length(est$piv)
	    out = complk(Y,R,rep(1,n),piv,Pi,Psi,k)
	    Phi = out$Phi; L = out$L; pv = out$pv
	   	V = array(0,c(n,k,TT))
   		Yvp = matrix(1/pv,n,k)
  		M = matrix(1,n,k)
   		V[,,TT] = Yvp*L[,,TT]
   		for(t in seq(TT-1,2,-1)){
   			M = (Phi[,,t+1]*M)%*%t(Pi[,,t+1])
      		V[,,t] = Yvp*L[,,t]*M
    	}
    	M = (Phi[,,2]*M)%*%t(Pi[,,2])
    	V[,,1] = Yvp*L[,,1]*M
# local deconding
   		Ul = matrix(0,n,TT)
		for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(V[i,,t])
    	if(n==1) Ul = as.vector(Ul)
# global deconding (Viterbi)
		R = L; Ug = matrix(0,n,TT)
		for(i in 1:n) for(t in 2:TT) for(u in 1:k) R[i,u,t] = Phi[i,u,t]*max(R[i,,t-1]*Pi[,u,t])
		if(n==1) Ug[,TT] = which.max(R[,,TT])
		else Ug[,TT] = apply(R[,,TT],1,which.max)
		for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*Pi[,Ug[i,t+1],t+1])
    	if(n==1) Ug = as.vector(Ug)
	}
	
# est_lm_cov_latent
	if(class(est)=="LMlatent"){
		param = est$param
		miss = any(is.na(Y))
		if(miss){
        	R = 1 * (!is.na(Y))
         	Y[is.na(Y)] = 0
		}else{
			R = NULL
		}
		if(dim(est$Psi)[3]==1){
			if(is.vector(Y)){
				Y = t(Y)
				if(miss) R = t(R)
				if(is.vector(X1)) X1 = t(X1) 		
				if(is.matrix(X2)) X2 = array(X2,c(1,dim(X2))) 
				if(is.vector(X2)) X2 = array(X2,c(1,length(X2),1))
			}		
			n = nrow(Y); TT = ncol(Y)
		}else{
			if(is.matrix(Y)){ 
				Y = array(Y,c(1,dim(Y)))
				if(miss) R = array(R,c(1,dim(R)))
				if(is.vector(X1)) X1 = t(X1)		
				if(is.matrix(X2)) X2 = array(X2,c(1,dim(X2)))	
				if(is.vector(X2)) X2 = array(X2,c(1,length(X2),1))
			}	
			n = dim(Y)[1]; TT = dim(Y)[2]; r = dim(Y)[3]
		}
		
			
		k = ncol(est$Be)+1
		Psi = est$Psi
		if(is.vector(X1)) X1 = matrix(X1,n,1)
		nc1 = dim(X1)[2] # number of covariates on the initial probabilities		
		Xlab = 1:n
		if(k == 2){
			GBe = as.matrix(c(0,1))
		}else{
			GBe = diag(k); GBe = GBe[,-1]
		}
		XXdis = array(0,c(k,(k-1)*(nc1+1),n))
		for(i in 1:n){
			xdis = c(1,X1[i,])
			XXdis[,,i] = GBe%*%(diag(k-1)%x%t(xdis))
		}
		be = as.vector(est$Be)
		out = prob_multilogit(XXdis,be,Xlab,fort)
		Piv = out$P
		if(is.matrix(X2)) X2 = array(X2,c(n,TT-1,1))
		nc2 = dim(X2)[3] # number of covariates on the transition probabilities
		Z = NULL
		for(t in 1:(TT-1)) Z = rbind(Z,X2[,t,])
		if(nc2==1) Z = as.matrix(X2)
		Zlab = 1:(n*(TT-1)); Zndis = n*(TT-1)
		if(param=="multilogit"){
   	 	ZZdis = array(0,c(k,(k-1)*(nc2+1),Zndis,k))
		    for(h in 1:k){
			    if(k==2){
				    if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
			    }else{
				    GGa = diag(k); GGa = GGa[,-h]
			    }  		
			    for(i in 1:Zndis){
				    zdis = c(1,Z[i,])
				    ZZdis[,,i,h] = GGa%*%(diag(k-1)%x%t(zdis))
			    }
		    }
		 }else if(param=="difflogit"){
   	     Zlab = (((Zlab-1)*k)%x%rep(1,k))+rep(1,n*(TT-1))%x%(1:k)
   	     ZZdis = array(0,c(k,k*(k-1)+(k-1)*nc2,Zndis*k))
   	     j = 0
			for(i in 1:Zndis){
   	         for(h in 1:k){
   	             j = j+1
   	             if(k==2){
				      if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
			        }else{
				        GGa = diag(k); GGa = GGa[,-h]
			        }  		
				    u = matrix(0,1,k); u[1,h] = 1
				    U = diag(k); U[,h] = U[,h]-1
				    U = U[,-1]
			        ZZdis[,,j] = cbind(u%x%GGa,U%x%t(Z[i,]))   
			       
   	         	}
		    }
		}
        if(param=="multilogit"){
    		Ga = matrix(est$Ga,(nc2+1)*(k-1),k)
	    	PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,n,TT))
		    for(h in 1:k){
			    out = prob_multilogit(ZZdis[,,,h],Ga[,h],Zlab,fort)
			    PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
		    }
		}else if(param=="difflogit"){
            Ga = c(as.vector(t(est$Ga[[1]])),as.vector(est$Ga[[2]]))
            PI = array(0,c(k,k,n,TT))
            out = prob_multilogit(ZZdis,Ga,Zlab,fort)
   		    PIdis = out$Pdis; 
   		    Tmp = array(out$P,c(k,n,TT-1,k))
	    	PI[,,,2:TT] = aperm(Tmp,c(1,4,2,3)) 
		}
		out = lk_comp_latent(Y,R,rep(1,n),Piv,PI,Psi,k,fort=fort)
		Phi = out$Phi; L = out$L; pv = out$pv

		out = prob_post_cov(Y,rep(1,n),Psi,Piv,PI,Phi,L,pv,fort=fort)
		V = out$V

# local deconding
   		Ul = matrix(0,n,TT)
		for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(V[i,,t])
    	if(n==1) Ul = as.vector(Ul)
# global deconding (Viterbi)
		R = L; Ug = matrix(0,n,TT)
		for(i in 1:n) for(t in 2:TT) for(u in 1:k) R[i,u,t] = Phi[i,u,t]*max(R[i,,t-1]*PI[,u,i,t])
		if(n==1) Ug[,TT] = which.max(R[,,TT])
		else Ug[,TT] = apply(R[,,TT],1,which.max)
		for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*PI[,Ug[i,t+1],i,t+1])
    	if(n==1) Ug = as.vector(Ug)		
	}

# est_lm_mixed
	if(class(est)=="LMmixed"){
		if(dim(est$Psi)[3]==1){
			if(is.vector(Y)) Y = t(Y)
			if(is.matrix(Y)) Y =  array(Y,c(dim(Y),1))	
			n = nrow(Y); TT = ncol(Y); r=1
		}else{
			if(is.matrix(Y)) Y = array(Y,c(1,dim(Y)))
			n = dim(Y)[1]; TT = dim(Y)[2]; r = dim(Y)[3]
		}
		yv = rep(1,n)
		la = est$la; Piv = est$Piv; Pi = est$Pi; Psi = est$Psi
		k1 = length(la); k2 = nrow(Piv)		
		Fc1 = matrix(0,n,k1); Fc2 = array(0,c(n,k1,k2)); Fc3 = array(0,c(n,k1,k2,k2))
		PP1 = array(0,c(n,k1,k2,TT))
		Phi = array(1,c(n,k2,TT))
		for(t in 1:TT){
  			#if(r==1) Phi[,,t] = Phi[,,t]*Psi[Y[,t]+1,,1] 
  			#else for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[Y[,t,j]+1,,j]
  			for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[Y[,t,j]+1,,j]
		}
		for(i in 1:n) for(u in 1:k1){
	   		 o = .Fortran("BWforback", TT, k2, Phi[i,,], Piv[,u], Pi[,,u], lk=0, Pp1=matrix(0,k2,TT), 
	        	         Pp2=array(0,c(k2,k2,TT)))
			Fc1[i,u] = exp(o$lk); Fc2[i,u,] = o$Pp1[,1]; Fc3[i,u,,] = apply(o$Pp2,c(1,2),sum)
  			PP1[i,u,,] = o$Pp1
  		}
		Fj1 = Fc1%*%diag(la)
		fm = rowSums(Fj1)
		fm = pmax(fm,10^-300)
		W = (Fj1/matrix(fm,n,k1))*yv
		U1 = apply(W,1,which.max)
		PV = array(W,c(n,k1,k2,TT))*PP1
# local decoding
   		Ul = matrix(0,n,TT)
		for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(PV[i,U1[i],,t])
    	if(n==1) Ul = as.vector(Ul)
# global deconding (Viterbi)
		R = array(0,c(n,k2,TT))
		for(i in 1:n) for(v in 1:k2) R[i,v,1] = Phi[i,v,1]*Piv[v,U1[i]]		 
		for(i in 1:n) for(t in 2:TT) for(v in 1:k2) R[i,v,t] = Phi[i,v,t]*max(R[i,,t-1]*Pi[,v,U1[i]])
		Ug = matrix(0,n,TT)
		if(n==1) Ug[,TT] = which.max(R[,,TT])
		else Ug[,TT] = apply(R[,,TT],1,which.max)
		for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*Pi[,Ug[i,t+1],U1[i]])
    	if(n==1) Ug = as.vector(Ug)
	}
	
#est_lm_manifest
if(class(est)=="LMmanifest"){
	lev = max(S)+1
	q = 1; k = dim(est$PI)[1]
	nt = prod(lev)
	n1=FALSE
	S = Y; X = X1
	if(is.vector(S)){
		n1 = TRUE
		n = 2; TT = length(S)
		S = rbind(t(S),0)
		if(is.matrix(X)){
			X1 = array(0,c(2,nrow(X),ncol(X)))
			X1[1,,] = X
			X = X1
		}else if (is.vector(X)){
			X1 = array(0,c(2,length(X),1))
			X1[1,,] = X
			X = X1
		}
	}else{
		n = nrow(S); TT = ncol(S)
	}

	if(is.array(S)) S = matrix(S,n,TT)
	if(is.matrix(X)) X = array(X,c(n,TT,1))
	
	Y0 = S+1
	S = array(0,c(nt,n,TT))
	for(i in 1:n) for(t in 1:TT){
   		ind = Y0[i,t]
   		S[ind,i,t] = 1
	}	

	nc = dim(X)[3]
	ne = lev-1

	XX = X
	X = array(0,c(ne,nc,n,TT))
	for(i in 1:n) for(t in 1:TT){
   		X[,,i,t] = rep(1,ne)%o%XX[i,t,]
	}

	out = marg_param(lev,"g")
	Cm = out$C; Mm = out$M
	Gm = cbind(-rep(1,lev-1),diag(lev-1))
	Hm = rbind(rep(0,lev-1),diag(lev-1))
	GHt = t(Gm)%*%t(Hm)
	lm = c(1,rep(0,lev-1))
	Lm = rbind(rep(0,lev-1),diag(lev-1))-rbind(diag(lev-1),rep(0,lev-1))

	if(q==1) sup = 0 else{
  		lim = 5;
  		sup = seq(-lim,lim,2*lim/(q-1))
	}
	Mar = diag(k)%x%matrix(1,1,q) 

 
	G2 = NULL; H2 = NULL; IPI = NULL
	if(k>1){
		for(c in 1:k){
    		G2c = diag(k)[,-c]   
      		H2c = diag(k)[-c,]; if (k==2) H2c[c]=-1 else H2c[,c]= -1
      		if(is.null(G2)) G2 = G2c else if(k==2) G2 = blkdiag(matrix(G2,ncol=1),matrix(G2c,ncol=1)) else G2 = blkdiag(G2,G2c) 
      		if(is.null(H2)) H2 = H2c else if(k==2) H2 = blkdiag(matrix(H2,nrow=1),matrix(H2c,nrow=1))else H2 = blkdiag(H2,H2c)
      		IPI = c(IPI,c+seq(0,k*(k-1),k))
    	}
	}

	mu = est$mu+est$al[1]
	al = est$al[-1]-est$al[1]
	la = est$la
	be = est$be
	PI = est$PI
	si = 3
	par = c(mu,al,si,be)
	if(k==1) tau = NULL else{
		tau = H2%*%log(PI[IPI]) 
	}
	las = la
	PIs = PI

	# find non-redundant X configurations (may be very slow)
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
	#Xd = Xd[,,1:nd]  #commentato per consentire 1 covariata
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
      			for(t in 1:TT) if(n==1) Pio[,j,t] = sum(S[,,t]*Pv1[,indn[,t]]) else Pio[,j,t] = colSums(S[,,t]*Pv1[,indn[,t]])
    		}
  		}
  		Q = rec1(Pio,las,PIs)
  
 		 if(q*k==1) pim = Q[,,TT] else if(n==1) pim = sum(Q[,,TT]) else pim = rowSums(Q[,,TT])
  		lk = sum(log(pim))
   
  		# E-step
  		out = rec3(Q,PIs,Pio,pim)
  		U = out$U; V = out$V
  		# M-step: latent parameters
  		if(k>1){
  			u1 = Mar%*%rowSums(U[,,1])
    		V1 = Mar%*%V%*%t(Mar)
    		out = lk_sta(tau,as.vector(u1),V1,G2,outl=TRUE)			
    		flk = out$flk; la = out$la; PI = out$PI
  		}
  		las = la; PIs = PI
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
    			for(h2 in 1:nd) XXRi1[,,h2] = cbind(D,Xd[,,h2])
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
        		XRi = (XXRi1[,,jd]%*%RRtc1[,,jd])%*%GHt
        		if(length(ind)==1){
        			s = s+XRi%*%(YGP[,ind]*pci)
        		}else{
            		s = s+XRi%*%(YGP[,ind]%*%pci)
        		}
        		FF = FF+sum(pci)*(XRi%*%Om[,,jd])%*%t(XRi)
     		}
   		}
  	}
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
  
  # Newton-Rapshon
  par1 = NULL;
  if(k>1) par1 = tau
  par1 = c(par1,par)
   
}
# separate parameters and compute aic and bic
mu = par[1:ne]
al = 0
if(k>1) al = c(al,par[(ne+1):(ne+k-1)])
mu = mu+al%*%la
al = al-al%*%la
be = par[(ne+k+1):length(par)]
np = k*(k-1)
np = np + (ne+(k-1)+nc) + ((k+1)*(q>1))
if(q==1){
	si=NULL; rho = NULL
}

# compute aic, bic and prediction of latent structure
out = lk_obs_manifest(par1,S,Xd,indn,lev,k,sup,G2,IPI,mod=0,outp=TRUE)
lk = out$lk; U = out$U
sup1 = t(Mar)%*%al
if(q>1) sup1 = sup1+matrix(1,k,1)%x%(sup*si)
PRED0 = array(0,c(n,k,TT)); PRED1 = matrix(0,n,TT)
for(t in 1:TT){
	PRED0[,,t] = U[,,t]%*%t(Mar)
    PRED1[,t] = U[,,t]%*%sup1
}

if(n1){
	V = array(PRED0[1,,],c(1,k,TT)); Phi = array(Pio[1,,],c(1,k,TT)) 
}else{ 
	V = PRED0; Phi = Pio
}
n = dim(V)[1]
piv = la
PI = PI
# local deconding
Ul = matrix(0,n,TT)
for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(V[i,,t])
if(n==1) Ul = as.vector(Ul)
# global deconding (Viterbi)
R = array(0,c(n,k,TT))
for(i in 1:n) for(v in 1:k) R[i,v,1] = Phi[i,v,1]*piv[v]	
Ug = matrix(0,n,TT)
for(i in 1:n) for(t in 2:TT) for(u in 1:k) R[i,u,t] = Phi[i,u,t]*max(R[i,,t-1]*PI[,u])
if(n==1) Ug[,TT] = which.max(R[,,TT])
else Ug[,TT] = apply(R[,,TT],1,which.max)
for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*PI[,Ug[i,t+1]])
if(n==1) Ug = as.vector(Ug)
	
}
# output
	out = list(Ul=Ul,Ug=Ug)
	
}
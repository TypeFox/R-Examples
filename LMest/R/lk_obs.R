lk_obs <- function(th,m,Bm,Cm,bv,k,S,R,yv,TT,r,mod){
	
# copute corresponding parameters
	miss=!is.null(R)
	b = max(bv)
    ns = dim(S)[1]
    n = sum(yv)
	Psi = array(NA,c(b+1,k,r))
	for(j in 1:r) Psi[1:(bv[j]+1),,j] = 0 
	for(u in 1:k) for(j in 1:r){
		indj = 1:(bv[j]+1)
		th1 = th[1:bv[j]]; th = th[-(1:bv[j])]
		Psi[indj,u,j] = exp(m[[j]]$Am%*%th1); Psi[indj,u,j] = Psi[indj,u,j]/sum(Psi[indj,u,j])				
	}
	th1 = th[1:(k-1)]; th = th[-(1:(k-1))]
	piv = exp(Bm%*%th1); piv = as.vector(piv/sum(piv))
	if(mod==0){
		Pi = array(0,c(k,k,TT))
		for(t in 2:TT) for(u in 1:k){
			th1 = th[1:(k-1)]; th = th[-(1:(k-1))]
			if(k==2) Pi[u,,t] = exp(Cm[,,u]*th1)
			else Pi[u,,t] = exp(Cm[,,u]%*%th1)
			Pi[u,,t] = Pi[u,,t]/sum(Pi[u,,t])
		}
	}
	if(mod==1){
		Pi = matrix(0,k,k)
		for(u in 1:k){
			th1 = th[1:(k-1)]; th = th[-(1:(k-1))]
			if(k==2) Pi[u,] = exp(Cm[,,u]*th1)
			else Pi[u,] = exp(Cm[,,u]%*%th1)
			Pi[u,] = Pi[u,]/sum(Pi[u,])
		}
		Pi = array(Pi,c(k,k,TT))
	}
	if(mod>1){
		Pi = array(0,c(k,k,TT))
		for(u in 1:k){
			th1 = th[1:(k-1)]; th = th[-(1:(k-1))]
			if(k==2) Pi[u,,2:mod] = exp(Cm[,,u]*th1)
			else Pi[u,,2:mod] = exp(Cm[,,u]%*%th1)
			Pi[u,,2:mod] = Pi[u,,2]/sum(Pi[u,,2])
		}
		for(u in 1:k){
			th1 = th[1:(k-1)]; th = th[-(1:(k-1))]
			if(k==2) Pi[u,,(mod+1):TT] = exp(Cm[,,u]*th1)
			else Pi[u,,(mod+1):TT] = exp(Cm[,,u]%*%th1)
			Pi[u,,(mod+1):TT] = Pi[u,,mod+1]/sum(Pi[u,,mod+1])
		}
	}
	Pi[,,1] = 0
# compute log-likelihood
    out = complk(S,R=R,yv,piv,Pi,Psi,k)
  	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
  	sc = NULL
# ---- E-step ----
# Compute V and U
   	V = array(0,c(ns,k,TT)); U = array(0,c(k,k,TT))
   	Yvp = matrix(yv/pv,ns,k)
  	M = matrix(1,ns,k);
   	V[,,TT] = Yvp*L[,,TT]
   	U[,,TT] = (t(L[,,TT-1])%*%(Yvp*Phi[,,TT]))*Pi[,,TT]
   	for(t in seq(TT-1,2,-1)){
   		M = (Phi[,,t+1]*M)%*%t(Pi[,,t+1]);
    	V[,,t] = Yvp*L[,,t]*M
    	U[,,t] = (t(L[,,t-1])%*%(Yvp*Phi[,,t]*M))*Pi[,,t]
    }
    M = (Phi[,,2]*M)%*%t(Pi[,,2])
    V[,,1] = Yvp*L[,,1]*M
# ---- M-step ----
# Update Psi
   	Y = array(0,c(b+1,TT,k,r))
   	for(j in 1:r) for(t in 1:TT) for(jb in 0:b){
    	ind = which(S[,t,j]==jb)
   		li = length(ind)
   		if(miss){
   			if(li==1) Y[jb+1,t,,j] = V[ind,,t]*R[ind,t,j]
			if(li>1) Y[jb+1,t,,j] = colSums(V[ind,,t]*R[ind,t,j])
   		}else{
    		if(li==1) Y[jb+1,t,,j] = V[ind,,t]
			if(li>1) Y[jb+1,t,,j] = colSums(V[ind,,t])
		}
	}
 	Y1 = apply(Y,c(1,3,4),sum)
 	#if(r==1) for(u in 1:k) sc = c(sc,t(Am)%*%(Y1[,u,1]-sum(Y1[,u,1])*Psi[,u,1]))
 	for(u in 1:k) for(j in 1:r){
 		indj = 1:(bv[j]+1)
 		sc = c(sc,t(m[[j]]$Am)%*%(Y1[indj,u,j]-sum(Y1[indj,u,j])*Psi[indj,u,j]))
 	}
# Update piv and Pi
	sc = c(sc,t(Bm)%*%(colSums(V[,,1])-n*piv))
	if(mod==0) {
	   	for(t in 2:TT) for(u in 1:k) sc = c(sc,t(Cm[,,u])%*%(U[u,,t]-sum(U[u,,t])*Pi[u,,t]))
	}
	if(mod==1){
	   	Ut = apply(U[,,2:TT],c(1,2),sum)
   	   	for(u in 1:k) sc = c(sc,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,2]))
	}
	if(mod>1){
		Ut1 = U[,,2:mod]
        if(length(dim(Ut1))>2) Ut1 = apply(Ut1,c(1,2),sum)
       	Ut2 = U[,,(mod+1):TT]
	    if(length(dim(Ut2))>2) Ut2 = apply(Ut2,c(1,2),sum)
    	for(u in 1:k) sc = c(sc,t(Cm[,,u])%*%(Ut1[u,]-sum(Ut1[u,])*Pi[u,,2]))
    	for(u in 1:k) sc = c(sc,t(Cm[,,u])%*%(Ut2[u,]-sum(Ut2[u,])*Pi[u,,mod+1]))
    	 
	}
# return
	out = list(lk=lk,sc=sc)	
}

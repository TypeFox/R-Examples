est_lm_basic <-
function(S,yv,k,start=0,mod=0,tol=10^-8,maxit=1000,out_se=FALSE,piv=NULL,Pi=NULL,Psi=NULL){

# Preliminaries
    check_der = FALSE  # to check derivatives
	n = sum(yv)
   	sS = dim(S)
  	ns = sS[1]
  	TT = sS[2]
  	if(min(S,na.rm=T)>0){
  		cat("|------------------- WARNING -------------------|\n")
  		cat("|The first response category must be coded as 0 |\n")
  		cat("|-----------------------------------------------|\n")	
 	}
 	 
  	if(is.data.frame(S)){
  		warning("Data frame not allowed for S")
  	}
  	
   if(ns!=length(yv)) stop("dimensions mismatch between S and yv")
  	
  	if(length(sS)==2){
  		r = 1
  		if(is.matrix(S)) S = array(S,c(dim(S),1))
  	}else r = sS[3]
  	
  	miss = any(is.na(S))
	if(miss){
         cat("Missing data in the dataset, treated as missing at random\n")
         R = 1 * (!is.na(S))
         S[is.na(S)] = 0
	}else{
		R = NULL
	}
  	Sv = matrix(S,ns*TT,r)
	if(miss) Rv = matrix(R,ns*TT,r)
 	bv = apply(Sv,2,max)
  	b = max(bv)
  	m = vector("list",r)
  	for(j in 1:r){
	    m[[j]]$Co = cbind(-diag(bv[j]),diag(bv[j]))
  		Maj = cbind(lower.tri(matrix(1,bv[j],bv[j]), diag = TRUE),rep(0,bv[j]))
  		m[[j]]$Ma = rbind(Maj,1-Maj)
  	}
  	th = NULL; sc = NULL; J = NULL  		
  	if(out_se){
  		for(j in 1:r){
	  		m[[j]]$A = cbind(-rep(1,bv[j]),diag(bv[j]))
  			m[[j]]$Am = rbind(rep(0,bv[j]),diag(bv[j]))
  		}
  		B = cbind(-rep(1,k-1),diag(k-1))
  		Bm = rbind(rep(0,k-1),diag(k-1))
  		C = array(0,c(k-1,k,k))
  		Cm = array(0,c(k,k-1,k))
  		for(u in 1:k){
  			C[,,u] = rbind(cbind(diag(u-1),-rep(1,u-1),matrix(0,u-1,k-u)),
  						   cbind(matrix(0,k-u,u-1),-rep(1,k-u),diag(k-u)))
  			Cm[,,u] = rbind(cbind(diag(u-1),matrix(0,u-1,k-u)),
  							rep(0,k-1),
	 					    cbind(matrix(0,k-u,u-1),diag(k-u)))
  		}
  	}
# When there is just 1 latent class
  	if(k == 1){
	    piv = 1; Pi = 1
   	 	P = matrix(NA,b+1,r)
   	 	for(j in 1:r) P[1:(bv[j]+1),j] = 0
	    for(t in 1:TT){
	      	for(j in 1:r){
		        	for(y in 0:b){
		        	ind = which(S[,t,j]==y)
		  		    P[y+1,j] = P[y+1,j]+sum(yv[ind])
		        }
	   		}
	    }
	    Psi = P/(n*TT)
	    pm = rep(1,ns)
	   	for(t in 1:TT) for(j in 1:r) pm = pm*Psi[S[,t,j]+1,j]
	    lk = sum(yv*log(pm))
	    np = r*b
	    aic = -2*lk+np*2
	    bic = -2*lk+np*log(n)
    		out = list(lk=lk,piv=piv,Pi=Pi,Psi=Psi,np=np,aic=aic,bic=bic,lkv=NULL,J=NULL,V=NULL,th=NULL,sc=NULL,call=match.call())
		class(out)="LMbasic" 
	    return(out)
  	}
# Starting values
	if(start == 0){
   		P = matrix(NA,b+1,r); E = matrix(NA,b,r)
   		for(j in 1:r) P[1:(bv[j]+1),j] = 0
        for(t in 1:TT) for(j in 1:r) for(y in 0:b){
	    		ind = which(S[,t,j]==y)
    			P[y+1,j] = P[y+1,j]+sum(yv[ind])
				E[1:bv[j],j] = m[[j]]$Co%*%log(m[[j]]$Ma%*%P[1:(bv[j]+1),j])
		}
  	   	Psi = array(NA,c(b+1,k,r)); Eta = array(NA,c(b,k,r))
        	grid = seq(-k,k,2*k/(k-1))
        	for(c in 1:k) for(j in 1:r){
			etac = E[1:bv[j],j]+grid[c]
			Eta[1:bv[j],c,j] = etac
			Psi[1:(bv[j]+1),c,j] = invglob(etac)
       	}
  		piv = rep(1,k)/k
		Pi = matrix(1,k,k)+9*diag(k); Pi = diag(1/rowSums(Pi))%*%Pi;
		Pi = array(Pi,c(k,k,TT)); Pi[,,1] = 0
  	}
  	if(start==1){
  		Psi = array(NA,c(b+1,k,r)) 
		for(j in 1:r){
			Psi[1:(bv[j]+1),,j] = matrix(runif((bv[j]+1)*k),bv[j]+1,k) 
			for(c in 1:k) Psi[1:(bv[j]+1),c,j] = Psi[1:(bv[j]+1),c,j]/sum(Psi[1:(bv[j]+1),c,j])
		}
	    Pi = array(runif(k^2*TT),c(k,k,TT))
	    for(t in 2:TT) Pi[,,t] = diag(1/rowSums(Pi[,,t]))%*%Pi[,,t]
	    Pi[,,1] = 0
	    piv = runif(k); piv = piv/sum(piv)
	}
	if(start==2){
		if(is.null(piv)) stop("initial value of the initial probabilities (piv) must be given in input")
		if(is.null(Pi)) stop("initial value of the transition probabilities (Pi) must be given in input")
		if(is.null(Psi)) stop("initial value of the conditional response probabilities (Psi) must be given in input")
		piv = piv
		Pi = Pi
		Psi = Psi 
	}
# Compute log-likelihood
    out = complk(S,R,yv,piv,Pi,Psi,k)
  	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
  	lk0 = sum(yv*log(yv/n)); dev = 2*(lk0-lk)
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat("     mod    |      k      |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
  	cat(sprintf("%11g",c(mod,k,start,0,lk)),"\n",sep=" | ")
  	it = 0; lko = lk-10^10; lkv = NULL
  	par = c(piv,as.vector(Pi),as.vector(Psi))
  	if(any(is.na(par))) par = par[-which(is.na(par))]
  	paro = par
# Iterate until convergence
	while((lk-lko)/abs(lk)>tol & it<maxit){
		Psi0 = Psi; piv0 = piv; Pi0 = Pi
		it = it+1;
# ---- E-step ----
# Compute V and U
#time = proc.time()
   		V = array(0,c(ns,k,TT)); U = array(0,c(k,k,TT))
   		Yvp = matrix(yv/pv,ns,k)
  		M = matrix(1,ns,k)
   		V[,,TT] = Yvp*L[,,TT]
   		U[,,TT] = (t(L[,,TT-1])%*%(Yvp*Phi[,,TT]))*Pi[,,TT]
   		for(t in seq(TT-1,2,-1)){
   			M = (Phi[,,t+1]*M)%*%t(Pi[,,t+1]);
      		V[,,t] = Yvp*L[,,t]*M
      		U[,,t] = (t(L[,,t-1])%*%(Yvp*Phi[,,t]*M))*Pi[,,t]
		}
		M = (Phi[,,2]*M)%*%t(Pi[,,2])
		V[,,1] = Yvp*L[,,1]*M
#print(proc.time()-time)
# If required store parameters
# ---- M-step ----
# Update Psi
	Y1 = array(NA,c(b+1,k,r))
	for(j in 1:r) Y1[1:(bv[j]+1)] = 0
	Vv = matrix(aperm(V,c(1,3,2)),ns*TT,k)
	for(j in 1:r) for(jb in 0:bv[j]) {
		ind = which(Sv[,j]==jb)
		if(length(ind)==1){
			if(miss) Y1[jb+1,,j] = Vv[ind,]*Rv[ind,j]
			else Y1[jb+1,,j] = Vv[ind,]
		}
		if(length(ind)>1){
			if(miss) Y1[jb+1,,j] = colSums(Vv[ind,]*Rv[ind,j])
			else Y1[jb+1,,j] = colSums(Vv[ind,])
		}
						
	}
	for(j in 1:r) for(c in 1:k){
		tmp = Y1[1:(bv[j]+1),c,j]
		tmp = pmax(tmp/sum(tmp),10^-10)
		Psi[1:(bv[j]+1),c,j] = tmp/sum(tmp)
	} 
#print(proc.time()-time)
# Update piv and Pi
	piv = colSums(V[,,1])/n
	U = pmax(U,10^-300)
	if(mod==0) for(t in 2:TT) Pi[,,t] = diag(1/rowSums(U[,,t]))%*%U[,,t]
	if(mod==1){
	    	Ut = apply(U[,,2:TT],c(1,2),sum)
  	    	Pi[,,2:TT] = array(diag(1/rowSums(Ut))%*%Ut,c(k,k,TT-1))
    	}
    	if(mod>1){
	    	Ut1 = U[,,2:mod]
        	if(length(dim(Ut1))>2) Ut1 = apply(Ut1,c(1,2),sum)
        	Ut2 = U[,,(mod+1):TT]
	    if(length(dim(Ut2))>2) Ut2 = apply(Ut2,c(1,2),sum)
    	   	Pi[,,2:mod] = array(diag(1/rowSums(Ut1,2))%*%Ut1,c(k,k,mod-1))         
        Pi[,,(mod+1):TT] = array(diag(1/rowSums(Ut2,2))%*%Ut2,c(k,k,TT-mod))         
	}
#print(proc.time()-time)
# Compute log-likelihood
    	paro = par; par = c(piv,as.vector(Pi),as.vector(Psi))
   	  	if(any(is.na(par))) par = par[-which(is.na(par))]
    	lko = lk
    	out = complk(S,R,yv,piv,Pi,Psi,k)
    	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
    	if(it/10 == round(it/10)) cat(sprintf("%11g",c(mod,k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
    	lkv = c(lkv,lk)
#print(proc.time()-time)
	}
# Compute information matrix if required	
	if(out_se){
		th = NULL
		for(u in 1:k) for(j in 1:r) th = c(th,m[[j]]$A%*%log(Psi[1:(bv[j]+1),u,j]))
		th = c(th,B%*%log(piv))
		if(mod==0) for(t in 2:TT) for(u in 1:k) th = c(th,C[,,u]%*%log(Pi[u,,t]))
		if(mod==1) for(u in 1:k) th = c(th,C[,,u]%*%log(Pi[u,,2]))
		if(mod>1) {
			for(u in 1:k) th = c(th,C[,,u]%*%log(Pi[u,,2]))
			for(u in 1:k) th = c(th,C[,,u]%*%log(Pi[u,,mod+1]))
		}		
		lth = length(th)	
		out = recursions(S,R,yv,Psi,piv,Pi,k,lth,m,Bm,Cm,bv,mod)
		F1 = out$F1; F2 = out$F2; F1d = out$F1d; F2d = out$F2d
	   	
	   	sc = NULL
	 	Y = array(NA,c(b+1,TT,k,r))
		for(j in 1:r) Y[1:(bv[j]+1),,,j] = 0
   		for(j in 1:r) for(t in 1:TT) for(jb in 0:bv[j]){
      		ind = which(S[,t,j]==jb)
      		if(length(ind)==1){
	      		if(miss) Y[jb+1,t,,j] = F1[,t,ind]*R[ind,t,j]*yv[ind]
	    	  		else Y[jb+1,t,,j] = F1[,t,ind]*yv[ind]
	    	  	}
      		if(length(ind)>1){
	      		if(miss) Y[jb+1,t,,j] = F1[,t,ind]%*%(R[ind,t,j]*yv[ind])
	    	  		else Y[jb+1,t,,j] = F1[,t,ind]%*%yv[ind]
	    	  	}
   		}
   		Y1 = apply(Y,c(1,3,4),sum)
   		for(u in 1:k) for(j in 1:r){
   			indj = 1:(bv[j]+1)
   			sc = c(sc,t(m[[j]]$Am)%*%(Y1[indj,u,j]-sum(Y1[indj,u,j])*Psi[indj,u,j]))
   		}
   		Y2 = Y1

   		for(u in 1:k) for(j in 1:r){
   			indj = 1:(bv[j]+1)
   			Juj = sum(Y1[indj,u,j])*t(m[[j]]$Am)%*%(diag(Psi[indj,u,j])-Psi[indj,u,j]%o%Psi[indj,u,j])%*%m[[j]]$Am
   			if(u==1 && j==1) J = Juj
  			else J = blkdiag(J,Juj)
   		}   	   	
   	   	bv1 = F1[,1,]%*%yv
		sc = c(sc,t(Bm)%*%(bv1-sum(bv1)*piv))
		J = blkdiag(J,n*t(Bm)%*%(diag(piv)-piv%o%piv)%*%Bm)
		if(mod==0){
			for(t in 2:TT){
				Ut = 0
				for(i in 1:ns) Ut = Ut+yv[i]*F2[,,t,i]
				for(u in 1:k){
   					sc = c(sc,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,t]))
   					J = blkdiag(J,sum(Ut[u,])*t(Cm[,,u])%*%(diag(Pi[u,,t])-Pi[u,,t]%o%Pi[u,,t])%*%Cm[,,u])
   				}
   			}
   		}
		if(mod==1){
			Ut = 0
			for(i in 1:ns) for(t in 2:TT) Ut = Ut+yv[i]*F2[,,t,i]
			for(u in 1:k){
   				sc = c(sc,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,2]))
   				J = blkdiag(J,sum(Ut[u,])*t(Cm[,,u])%*%(diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2])%*%Cm[,,u])
   			}
   		}
   		if(mod>1){
   			Ut=0
   			for(i in 1:ns) for(t in 2:mod) Ut = Ut+yv[i]*F2[,,t,i] 
   			for(u in 1:k){
   				sc = c(sc,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,2]))
   				J = blkdiag(J,sum(Ut[u,])*t(Cm[,,u])%*%(diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2])%*%Cm[,,u])
   			}
   			Ut=0
   			for(i in 1:ns) for(t in (mod+1):TT) Ut = Ut+yv[i]*F2[,,t,i] 
   			for(u in 1:k){
   				sc = c(sc,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,mod+1]))
   				J = blkdiag(J,sum(Ut[u,])*t(Cm[,,u])%*%(diag(Pi[u,,mod+1])-Pi[u,,mod+1]%o%Pi[u,,mod+1])%*%Cm[,,u])
   			}
   			
   		}
   		J = as.matrix(J) 
   	    Jd = NULL
		for(pa in 1:lth){
			scj = NULL
	   		Y = array(NA,c(b+1,TT,k,r))
	   		for(j in 1:r) Y[1:(bv[j]+1),,,j] = 0
   			for(j in 1:r) for(t in 1:TT) for(jb in 0:b){
   		   		ind = which(S[,t,j]==jb)
   		   		if(length(ind)==1){
	   		   		if(miss) Y[jb+1,t,,j] = F1d[,t,ind,pa]*R[ind,t,j]*yv[ind]
		   	   		else Y[jb+1,t,,j] = F1d[,t,ind,pa]*yv[ind]
		   	   	}
		   	   	if(length(ind)>1){
	   		   		if(miss) Y[jb+1,t,,j] = F1d[,t,ind,pa]%*%(R[ind,t,j]*yv[ind])
		   	   		else Y[jb+1,t,,j] = F1d[,t,ind,pa]%*%yv[ind]
		   	   	}
   			}
	   		Y1 = apply(Y,c(1,3,4),sum)
	   	   	for(u in 1:k) for(j in 1:r){
	   	   		indj = 1:(bv[j]+1)
	   	   		scj = c(scj,t(m[[j]]$Am)%*%(Y1[indj,u,j]-sum(Y1[indj,u,j])*Psi[indj,u,j]))
	   	   	}
	   	   	bv1 = F1d[,1,,pa]%*%yv
			scj = c(scj,t(Bm)%*%(bv1-sum(bv1)*piv))
			if(mod==0) {
				for(t in 2:TT){
					Ut = 0
					for(i in 1:ns) Ut = Ut+yv[i]*F2d[,,t,i,pa]
				 	for(u in 1:k) scj = c(scj,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,t]))
				 }
			}
			if(mod==1){
				Ut = 0
				for(i in 1:ns) for(t in 2:TT) Ut = Ut+yv[i]*F2d[,,t,i,pa]
				for(u in 1:k) scj = c(scj,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,2]))
			}
			if(mod>1){
				Ut = 0
				for(i in 1:ns) for(t in 2:mod) Ut = Ut+yv[i]*F2d[,,t,i,pa]
				for(u in 1:k) scj = c(scj,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,2]))
				Ut = 0
				for(i in 1:ns) for(t in (mod+1):TT) Ut = Ut+yv[i]*F2d[,,t,i,pa]
				for(u in 1:k) scj = c(scj,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,mod+1]))
			}
			
			Jd = cbind(Jd,scj)
		}
		J = J-Jd
		Va = ginv(J)
		for(u in 1:k) for(j in 1:r){
			indj = 1:(bv[j]+1)
			Om = diag(Psi[indj,u,j])-Psi[indj,u,j]%o%Psi[indj,u,j]
			if(u==1) {
				if(j==1) M = Om%*%m[[j]]$Am
				else M = blkdiag(M,Om%*%m[[j]]$Am)
			} else M = blkdiag(M,Om%*%m[[j]]$Am)
		}
		Om = diag(piv)-tcrossprod(piv,piv)
		M = blkdiag(M,Om%*%Bm)
		if(mod==0){
			for(t in 2:TT) for(u in 1:k){
				Om = diag(Pi[u,,t])-Pi[u,,t]%o%Pi[u,,t]
				M = blkdiag(M,Om%*%Cm[,,u])
			}
		}
		if(mod==1){
			for(u in 1:k){
				Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
				M = blkdiag(M,Om%*%Cm[,,u])
			}
		}
		if(mod>1){
			for(u in 1:k){
				Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
				M = blkdiag(M,Om%*%Cm[,,u])
			}
			for(u in 1:k){
				Om = diag(Pi[u,,mod+1])-Pi[u,,mod+1]%o%Pi[u,,mod+1]
				M = blkdiag(M,Om%*%Cm[,,u])
			}
		}
		M = as.matrix(M)
		Va = M%*%Va%*%t(M)
		dVa = diag(Va)
		if(any(dVa<0)) warning("Negative elements in the estimated variance-covariance matrix for the parameters estimates")
		se = sqrt(abs(dVa))
		# Divide parameters
		nPsi = sum(bv+1)*k
		sePsi = se[1:nPsi]
		sepiv = se[nPsi+(1:k)]
		if(mod==0) sePi = se[nPsi+k+(1:(k*k*(TT-1)))]
		if(mod==1) sePi = se[nPsi+k+(1:(k*k))]
		if(mod>1) sePi = se[nPsi+k+(1:(k*k*2))]
#to check derivatives
		if(check_der){
			J0 = J
			th0 = th-10^-5/2
			out = lk_obs(th0,m,Bm,Cm,bv,k,S,R=R,yv,TT,r,mod)
			lk0 = out$lk; sc0 = out$sc
			lth = length(th)
			scn = rep(0,lth)
			J = matrix(0,lth,lth)
			for(j in 1:lth){
				thj = th0; thj[j] = thj[j]+10^-5
				out = lk_obs(thj,m,Bm,Cm,bv,k,S,R=R,yv,TT,r,mod)
				scn[j] = (out$lk-lk0)/10^-5
				J[,j] = (out$sc-sc0)/10^-5
			}
			J = -(J+t(J))/2
			print(c(lk,lk0))
			print(round(cbind(sc,scn,sc0),5))
			print(round(cbind(diag(J),diag(J0),diag(J)-diag(J0)),4))
			browser()
		}
	}
# Compute number of parameters  
	np = (k-1)+k*sum(bv)
  	if(mod==0) np = np+(TT-1)*k*(k-1)
  	if(mod==1) np = np+k*(k-1)
  	if(mod==2) np = np+2*k*(k-1)
  	aic = -2*lk+np*2
  	bic = -2*lk+np*log(n)
	cat(sprintf("%11g",c(mod,k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")	
	# adjust output
	lk = as.vector(lk)
	dimnames(Pi)=list(state=1:k,state=1:k,time=1:TT)
	dimnames(Psi)=list(category=0:b,state=1:k,item=1:r)

	out = list(lk=lk,piv=piv,Pi=Pi,Psi=Psi,np=np,aic=aic,bic=bic,lkv=lkv,V=V,call=match.call())
	if(out_se){
		sePsi0 = sePsi
		sePsi = array(NA,c(b+1,k,r))
		ind = 0
		for(u in 1:k) for(j in 1:r){
			indj = 1:(bv[j]+1)
			ind = max(ind)+indj
			sePsi[indj,u,j] = sePsi0[ind]
		}
		sePsi = array(sePsi,c(b+1,k,r))
		sePi0 = sePi
		sePi = array(0,c(k,k,TT))
		if(mod>1){
			sePi0 = array(sePi0,c(k,k,2))
			sePi0 = aperm(sePi0,c(2,1,3))
			sePi[,,2:mod] = sePi0[,,1]
			sePi[,,(mod+1):TT] = sePi0[,,2]
		} else {
			sePi[,,2:TT] = sePi0
			sePi = aperm(sePi,c(2,1,3))
		}
		dimnames(sePsi) = list(category=0:b,state=1:k,item=1:r)
		dimnames(sePi) = list(state=1:k,state=1:k,time=1:TT)
				
		out$sepiv = sepiv
		out$sePi = sePi
		out$sePsi = sePsi
	}
	
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    class(out)="LMbasic"
	return(out)
}

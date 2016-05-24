est_lm_mixed <- function(S,yv=rep(1,nrow(S)),k1,k2,start=0,tol=10^-8,maxit=1000,out_se=FALSE){

# EM algorithm for mixed HMM based
# Y = matrix of data
# k1 = number of latent classes
# k2 = number of latent states
# start = 0 for deterministic initialization, 1 for stochastic initialization
# tol = tollerance level to stop

# preliminaries
	k2 = as.integer(k2)
	n = sum(yv)
   	sS = dim(S)
  	ns = sS[1]
  	TT = sS[2]
  	if(is.data.frame(S)) warning("Data frame not allowed for S")
  
  	if(min(S)>0){
  		cat("|------------------------------------------ WARNING -----------------------------------------|\n")
  		cat("|The first response category must be coded as 0                                              |\n")
  		cat("|------------------------------------------ WARNING -----------------------------------------|\n")
 	} 
  	
    if(ns!=length(yv)) stop("dimensions mismatch between S and yv")    
  	if(length(sS)==2) r = 1
  	else r = sS[3]
  	if(r==1) {
  		if(is.matrix(S)) S = array(S,c(dim(S),1))
  	}	
	Sv = matrix(S,ns*TT,r)
	
  	b = max(S)
    flag = FALSE
    for (j in 1:r) if (length(table(Sv[, j])) < b) flag = TRUE
    if(flag){
   		cat("|------------------------------------------ WARNING -----------------------------------------|\n")
        cat("| Response variables must have the same number of categories                                 |\n")
		cat("|------------------------------------------ WARNING -----------------------------------------|\n")
    }
    
    Co = cbind(-diag(b),diag(b))
  	Ma = cbind(lower.tri(matrix(1,b,b), diag = TRUE),rep(0,b))
  	Ma = rbind(Ma,1-Ma)
# starting parameters
	if(start==0){
		la = rep(1,k1)/k1
   		Pi = matrix(1,k2,k2)+9*diag(k2); Pi = Pi/rowSums(Pi)
		Pi = array(Pi,c(k2,k2,k1))
   		P = matrix(0,b+1,r)
       	for(t in 1:TT) for(j in 1:r) for(y in 0:b){
			#if(r==1) ind = which(S[,t]==y) else ind = which(S[,t,j]==y)
			 ind = which(S[,t,j]==y)
    			P[y+1,j] = P[y+1,j]+sum(yv[ind])
    	}
		E = Co%*%log(Ma%*%P)
  	   	Psi = array(0,c(b+1,k2,r)); Eta = array(0,c(b,k2,r))
       	if(k2 == 1) grid = k2 else grid = seq(-k2,k2,2*k2/(k2-1))
       	for(c in 1:k2) for(j in 1:r){
     		etac = E[,j]+grid[c]
       		Eta[,c,j] = etac
       		Psi[,c,j] = invglob(etac)
       	}
	}else{
	    Psi = array(runif((b+1)*k2*r),c(b+1,k2,r)) 
	    for(j in 1:r) for(c in 1:k2) Psi[,c,j] = Psi[,c,j]/sum(Psi[,c,j])
		la = runif(k1); la = la/sum(la)
		Pi = array(0,c(k2,k2,k1))
		for(u in 1:k1){
			Pi[,,u] = matrix(runif(k2^2),k2,k2)
			Pi[,,u] = Pi[,,u]/rowSums(matrix(Pi[,,u],k2,k2))
		}
	}
	if(k1==1){
		if(start==0){
			Piv = matrix(1,k2,1)/k2
		}else{
			temp = runif(k2)
			Piv = matrix(temp/sum(temp),k2,1)
		}
	}else{
		if(start==0){
			if(k2==1){
				Piv = matrix(1,1,k1)				
			}else{
				Piv = matrix(0,k2,k1)	
				gl = log((k2-1):1)-log(1:(k2-1))
				for(u in 1:k1) Piv[,u] = diff(c(0,1/(1+exp(gl+u-(1+k1)/2)),1))
			}
		}else{
			Piv = matrix(0,k2,k1)	
			for(u in 1:k1){
				temp = runif(k2)
				Piv[,u] = temp/sum(temp)
			}
		}
	}

# EM algorithm for composite likelihood one-wise case row by row
# *** compute log-likelihood ***
# conditional distribution of any row given the row effect and corresponding joint
	Fc1 = matrix(0,ns,k1); Fc2 = array(0,c(ns,k1,k2)); Fc3 = array(0,c(ns,k1,k2,k2))
	PP1 = array(0,c(ns,k1,k2,TT))
	Phi = array(1,c(ns,k2,TT))
	for(t in 1:TT){
  		#if(r==1) Phi[,,t] = Phi[,,t]*Psi[S[,t]+1,,1] else for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[S[,t,j]+1,,j]
  		for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[S[,t,j]+1,,j]
	}
	for(i in 1:ns) for(u in 1:k1){
	    o = .Fortran("BWforback", TT, k2, Phi[i,,], Piv[,u], Pi[,,u], lk=0, Pp1=matrix(0,k2,TT), 
	                 Pp2=array(0,c(k2,k2,TT)))
		Fc1[i,u] = exp(o$lk); Fc2[i,u,] = o$Pp1[,1]; Fc3[i,u,,] = apply(o$Pp2,c(1,2),sum)
  		PP1[i,u,,] = o$Pp1
  	}
	Fj1 = Fc1%*%diag(la)
	fm = rowSums(Fj1)
	fm = pmax(fm,10^-300)
	lk = yv%*%log(fm)
	W = (Fj1/matrix(fm,ns,k1))*yv

# iterate until convergence
	it = 0; lko = -Inf

	cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat("     k1     |      k2     |    start    |     step    |     lk      |    lk-lko   |\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
  	cat(sprintf("%11g",c(k1,k2,start,0,lk)),"\n",sep=" | ")
 	
	while((lk-lko)/abs(lk)>tol & it<maxit){
		it = it+1
  		lko = lk

# # E-step update row-weights
		WZ1 = array(W,c(ns,k1,k2))*Fc2
		WZ2 = array(W,c(ns,k1,k2,k2))*Fc3
		PV = array(W,c(ns,k1,k2,TT))*PP1
		PV = apply(PV,c(1,3,4),sum)

# # M-step
		la = colSums(W); la = la/sum(la)
		for(u in 1:k1){
			Piv[,u] = apply(matrix(WZ1[,u,],ns,k2),2,sum)
			Piv[,u] = pmax(Piv[,u],10^-20)
			Piv[,u] = Piv[,u]/sum(Piv[,u])
			Pi[,,u] = apply(array(WZ2[,u,,],c(ns,k2,k2)),c(2,3),sum)
			Pi = pmax(Pi,10^-20)
			Pi[,,u] = Pi[,,u]/rowSums(matrix(Pi[,,u],k2,k2))
		}
		
		Y1 = array(0,c(b+1,k2,r))
		Vv = matrix(aperm(PV,c(1,3,2)),ns*TT,k2)
		for(j in 1:r) for(y in 0:b) {
			ind = which(Sv[,j]==y)
			if(k2==1) Y1[y+1,,j] = sum(Vv[ind,]) else Y1[y+1,,j] = colSums(Vv[ind,])				
		}
		for(j in 1:r) for(c in 1:k2) Psi[,c,j] = Y1[,c,j]/sum(Y1[,c,j]) 


# # *** compute log-likelihood ***
# # conditional distribution of any row given the row effect and corresponding joint
		Phi = array(1,c(ns,k2,TT))
		for(t in 1:TT){
	  		#if(r==1) Phi[,,t] = Phi[,,t]*Psi[S[,t]+1,,1]
	  		#else for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[S[,t,j]+1,,j]
	  		for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[S[,t,j]+1,,j]
		}
		for(i in 1:ns) for(u in 1:k1){
		    o = .Fortran("BWforback", TT, k2, Phi[i,,], Piv[,u], Pi[,,u], lk=0, Pp1=matrix(0,k2,TT), 
		                 Pp2=array(0,c(k2,k2,TT)))
			Fc1[i,u] = exp(o$lk); Fc2[i,u,] = o$Pp1[,1]; Fc3[i,u,,] = apply(o$Pp2,c(1,2),sum)
	  		PP1[i,u,,] = o$Pp1
	  	}
		Fj1 = Fc1%*%diag(la)
		fm = rowSums(Fj1)
		fm = pmax(fm,10^-300)
		lk = yv%*%log(fm)
		W = (Fj1/matrix(fm,ns,k1))*yv
# display output	
		if(it/10 == floor(it/10)) cat(sprintf("%11g",c(k1,k2,start,it,lk,lk-lko)),"\n",sep=" | ")
	}
	if(it/10 > floor(it/10))  cat(sprintf("%11g",c(k1,k2,start,it,lk,lk-lko)),"\n",sep=" | ")

# BIC
	np = (k1-1)+k1*(k2^2-1)+k2*r*b
	bic = -2*lk+np*log(n)
	
# compute standard errors if required
	if(out_se){
# structure of parameter vector
		lla = logit1(la)$lp
		Der = expit1(lla)$Der
		lPiv = matrix(0,k2-1,k1)
		for(u in 1:k1){
			lPiv[,u] = logit1(Piv[,u])$lp
			Der = blkdiag(Der,expit1(lPiv[,u])$Der)
		}
		lPiv = as.vector(lPiv)
		lPi = array(0,c(k2-1,k2,k1))
		for(u in 1:k1) for(v in 1:k2){
			lPi[,v,u] = logit1(Pi[v,,u],v)$lp
			Der = blkdiag(Der,expit1(lPi[,v,u])$Der)
		}
		lPi = as.vector(lPi)
		lPsi = array(0,c(b,k2,r))
		for(j in 1:r) for(v in 1:k2){
			lPsi[,v,j] = logit1(Psi[,v,j])$lp
			Der = blkdiag(Der,expit1(lPsi[,v,j])$Der)
		} 
		lPsi = as.vector(lPsi)
		th = c(lla,lPiv,lPi,lPsi)
		nla = length(lla); nPiv = length(lPiv); nPi = length(lPi); nPsi = length(lPsi)
		out = lk_obs_mixed(th,nla,nPiv,nPi,nPsi,S,yv,r,k1,k2)
		scn = NULL; Jn = NULL
		for(j in 1:length(th)){
			th1 = th; th1[j] = th1[j]+10^-6
			out1 = lk_obs_mixed(th1,nla,nPiv,nPi,nPsi,S,yv,r,k1,k2)	
			scn = c(scn,(out1$lk-out$lk)*10^6)	
			Jn = cbind(Jn,(out1$sc-out$sc)*10^6)	
		}
		Jn = (Jn+t(Jn))/2
		Vn = ginv(-Jn)
		if(any(diag(Vn)<0)) print("negative elements in the diagonal of inv(Information)")
		Vn1 = Der%*%Vn%*%t(Der)
		se1 = sqrt(abs(diag(Vn1)))
		sela = se1[1:k1]; se1 = se1[-(1:k1)]
		sePiv = matrix(0,k2,k1)
		for(u in 1:k1){
			sePiv[,u] = se1[1:k2]
			se1 = se1[-(1:k2)]
		}
		dimnames(sePiv)=list(v=1:k2,u=1:k1)
		sePi = array(0,c(k2,k2,k1))
		for(u in 1:k1) for(v in 1:k2){
			sePi[v,,u] = se1[1:k2]
			se1 = se1[-(1:k2)]
		}
		dimnames(sePi)=list(v0=1:k2,v1=1:k2,u=1:k1)
		sePsi = array(0,c(b+1,k2,r))
		for(j in 1:r) for(v in 1:k2){
			sePsi[,v,j] = se1[1:(b+1)]
			se1 = se1[-(1:(b+1))]
		} 
		dimnames(sePsi)=list(y=0:b,v=1:k2,j=1:r)
		# print(c(lk,out$lk,lk-out$lk))
		# print(cbind(out$sc,scn,out$sc-scn))
	}

# adjust output
    lk = as.vector(lk); bic = as.vector(bic)
	dimnames(Piv)=list(v=1:k2,u=1:k1)
	dimnames(Pi)=list(v0=1:k2,v1=1:k2,u=1:k1)
	dimnames(Psi)=list(y=0:b,v=1:k2,j=1:r)
	out = list(la=la,Piv=Piv,Pi=Pi,Psi=Psi,lk=lk,W=W,np=np,bic=bic,call=match.call())
	if(out_se){out$sela = sela; out$sePiv = sePiv; out$sePi = sePi; out$sePsi = sePsi}
	cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
	class(out)="LMmixed"
	out
	
}
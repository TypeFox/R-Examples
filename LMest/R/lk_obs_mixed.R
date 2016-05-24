lk_obs_mixed <- function(th,nla,nPiv,nPi,nPsi,S,yv,r,k1,k2){

# preliminaries
	b = max(S)
	ns = nrow(S)
	TT = ncol(S)
	Sv = matrix(S,ns*TT,r)	
# separate parameters
	lla = th[1:nla]; th = th[-(1:nla)]
	lPiv = th[1:nPiv]; th = th[-(1:nPiv)]
	lPi = th[1:nPi]; th = th[-(1:nPi)]
	lPsi = th[1:nPsi]; th = th[-(1:nPsi)]
	la = expit1(lla)$p		
	lPiv = matrix(lPiv,k2-1,k1)
	Piv = matrix(0,k2,k1)
	for(u in 1:k1) Piv[,u] = expit1(lPiv[,u])$p
	lPi = array(lPi,c(k2-1,k2,k1))
	Pi = array(0,c(k2,k2,k1))
	for(u in 1:k1) for(v in 1:k2) Pi[v,,u] = expit1(lPi[,v,u],v)$p
	lPsi = array(lPsi,c(b,k2,r))
	Psi = array(0,c(b+1,k2,r))
	for(j in 1:r) for(v in 1:k2) Psi[,v,j] = expit1(lPsi[,v,j])$p 
# compute log-likelihood
	Fc1 = matrix(0,ns,k1); Fc2 = array(0,c(ns,k1,k2)); Fc3 = array(0,c(ns,k1,k2,k2))
	PP1 = array(0,c(ns,k1,k2,TT))
	Phi = array(1,c(ns,k2,TT))
	for(t in 1:TT){
  		if(r==1) Phi[,,t] = Phi[,,t]*Psi[S[,t]+1,,1]
  		else for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[S[,t,j]+1,,j]
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
	lk = sum(yv*log(fm))	
	W = (Fj1/matrix(fm,ns,k1))*yv
# compute score
# E-step update row-weights
	WZ1 = array(W,c(ns,k1,k2))*Fc2
	WZ2 = array(W,c(ns,k1,k2,k2))*Fc3
	PV = array(W,c(ns,k1,k2,TT))*PP1
	PV = apply(PV,c(1,3,4),sum)
# M-step
	n = sum(yv)
	scla = (colSums(W)-n*la)[-1]
	scPiv = matrix(0,k2-1,k1)
	scPi = array(0,c(k2-1,k2,k1))
	for(u in 1:k1){
		tmp = apply(matrix(WZ1[,u,],ns,k2),2,sum)
		scPiv[,u] = (tmp-sum(tmp)*Piv[,u])[-1]
		Tmp = apply(array(WZ2[,u,,],c(ns,k2,k2)),c(2,3),sum)
		for(v in 1:k2) scPi[,v,u] = (Tmp[v,]-sum(Tmp[v,])*Pi[v,,u])[-v]
	}
	scPsi = array(0,c(b,v,r))
	Y1 = array(0,c(b+1,k2,r))
	Vv = matrix(aperm(PV,c(1,3,2)),ns*TT,k2)
	for(j in 1:r) for(y in 0:b) {
		ind = which(Sv[,j]==y)
		Y1[y+1,,j] = colSums(Vv[ind,])				
	}
	for(j in 1:r) for(c in 1:k2) scPsi[,c,j] = (Y1[,c,j]-sum(Y1[,c,j])*Psi[,c,j])[-1] 

# output
	sc = c(scla,scPiv,scPi,scPsi)
	out = list(lk=lk,sc=sc)
}
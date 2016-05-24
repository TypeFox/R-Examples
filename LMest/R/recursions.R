recursions <- function(S,R,yv,Psi,piv,Pi,k,lth,matr,Bm,Cm,bv,mod){
	
	
# forward recursion
sS = dim(S)
ns = sS[1]
TT = sS[2]
if(length(sS)==2){ 
	r = 1
	S = array(S,c(ns,TT,r))
}else r = sS[3]
b = max(bv)
miss = !is.null(R)
M = array(1,c(k,TT,ns))
if(miss) for(i in 1:ns) for(t in 1:TT) for(j in 1:r) M[,t,i] = M[,t,i]*(Psi[S[i,t,j]+1,,j]*R[i,t,j]+(1-R[i,t,j]))
else for(i in 1:ns) for(t in 1:TT) for(j in 1:r) M[,t,i] = M[,t,i]*Psi[S[i,t,j]+1,,j]

Q = array(0,c(k,TT,ns))
pv = rep(0,ns)
for(i in 1:ns){
	q = M[,1,i]*piv
	Q[,1,i] = q
	for(t in 2:TT){
		q = M[,t,i]*(t(Pi[,,t])%*%q)
		Q[,t,i] = q	
	}
	pv[i] = sum(q)	
}

#log-likelihood
lk = sum(yv*log(pv))

# backward recurion
Qb = array(0,c(k,TT,ns))
for(i in 1:ns){
	qb = rep(1,k)
	Qb[,TT,i] = qb
	for(t in seq(TT-1,1,-1)){
		qb = (Pi[,,t+1])%*%(M[,t+1,i]*qb)
		Qb[,t,i] = qb	
	}
}

# posterior probabilities
F1 = Q*Qb
for(i in 1:ns) F1[,,i] = F1[,,i]/pv[i]
F2 = array(0,c(k,k,TT,ns))
QbM = Qb*M
for(i in 1:ns) for(t in 2:TT) F2[,,t,i] = Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]


#derivative of the forward recursion
ind = 0
#derivative of Psi
Psid = array(NA,c(b+1,k,r,lth))
for(j in 1:r) Psid[1:(bv[j]+1),,j,] = 0
for(u in 1:k) for(j in 1:r){
	indj = 1:(bv[j]+1)
	Om = diag(Psi[indj,u,j])-Psi[indj,u,j]%o%Psi[indj,u,j]
	D = Om%*%matr[[j]]$Am
	for(h in 1:bv[j]){
		ind = ind+1
		Psid[indj,u,j,ind] = D[,h]
	}		
}
#derivative of piv
pivd = matrix(0,k,lth)
Om = diag(piv)-piv%o%piv
D = Om%*%Bm
for(j in 1:(k-1)){
	ind = ind+1
	pivd[,ind] = D[,j]
}
#derivative of Pi
if(mod==0){
	Pid = array(0,c(k,k,TT,lth))
	for(t in 2:TT) for(u in 1:k){
		Om = diag(Pi[u,,t])-Pi[u,,t]%o%Pi[u,,t]
		D = Om%*%Cm[,,u]
		for(j in 1:(k-1)){
			ind = ind+1
			Pid[u,,t,ind] = D[,j]
		}	
	}
}
if(mod==1){
	Pid = array(0,c(k,k,lth))
	for(u in 1:k){
		Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
		D = Om%*%Cm[,,u]
		for(j in 1:(k-1)){
			ind = ind+1
			Pid[u,,ind] = D[,j]
		}	
	}
}
if(mod>1){
	Pid = array(0,c(k,k,2,lth))
	for(u in 1:k){
		Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
		D = Om%*%Cm[,,u]
		for(j in 1:(k-1)){
			ind = ind+1
			Pid[u,,1,ind] = D[,j]
		}	
	}	
	for(u in 1:k){
		Om = diag(Pi[u,,mod+1])-Pi[u,,mod+1]%o%Pi[u,,mod+1]
		D = Om%*%Cm[,,u]
		for(j in 1:(k-1)){
			ind = ind+1
			Pid[u,,2,ind] = D[,j]
		}	
	}	
}

#derivative of the forward recursion
Md = array(0,c(k,TT,ns,lth))
for(i in 1:ns) for(t in 1:TT) for(j in 1:r){
	if(miss) Md[,t,i,] = Md[,t,i,]+(Psid[S[i,t,j]+1,,j,]*R[i,t,j])*M[,t,i]/(Psi[S[i,t,j]+1,,j]*R[i,t,j]+(1-R[i,t,j]))
	else Md[,t,i,] = Md[,t,i,]+Psid[S[i,t,j]+1,,j,]*M[,t,i]/Psi[S[i,t,j]+1,,j]
} 

Qd = array(0,c(k,TT,ns,lth))
pvd = matrix(0,ns,lth)
for(i in 1:ns){
	Qd[,1,i,] = Md[,1,i,]*piv+M[,1,i]*pivd
	for(t in 2:TT){
		q = Q[,t-1,i]
		if(mod==0) Qd[,t,i,] = Md[,t,i,]*as.vector(t(Pi[,,t])%*%q)+M[,t,i]*(t(Pi[,,t])%*%Qd[,t-1,i,])
		for(j in 1:lth){
			qd = Qd[,t-1,i,j]
			if(mod==0) Qd[,t,i,j] = Qd[,t,i,j]+M[,t,i]*(t(Pid[,,t,j])%*%q)
			if(mod==1) Qd[,t,i,j] = Md[,t,i,j]*(t(Pi[,,t])%*%q)+M[,t,i]*(t(Pid[,,j])%*%q)+M[,t,i]*(t(Pi[,,t])%*%qd)
			if(mod>1){
				if(t<=mod) Qd[,t,i,j] = Md[,t,i,j]*(t(Pi[,,t])%*%q)+M[,t,i]*(t(Pid[,,1,j])%*%q)+M[,t,i]*(t(Pi[,,t])%*%qd)
				if(t>mod) Qd[,t,i,j] = Md[,t,i,j]*(t(Pi[,,t])%*%q)+M[,t,i]*(t(Pid[,,2,j])%*%q)+M[,t,i]*(t(Pi[,,t])%*%qd)
			}
		}
	}
	for(j in 1:lth) pvd[i,j] = sum(Qd[,TT,i,j])	
}


#score
sc = (yv/pv)%*%pvd

# derivative of backward recurion
Qbd = array(0,c(k,TT,ns,lth))
for(j in 1:lth){
	for(i in 1:ns){
		qbd = rep(0,k)
		Qbd[,TT,i,j] = qbd
		for(t in seq(TT-1,1,-1)){
			qb = Qb[,t+1,i]
			if(mod==0) qbd = (Pid[,,t+1,j])%*%(M[,t+1,i]*qb)+(Pi[,,t+1])%*%(Md[,t+1,i,j]*qb)+(Pi[,,t+1])%*%(M[,t+1,i]*qbd)
			if(mod==1) qbd = (Pid[,,j])%*%(M[,t+1,i]*qb)+(Pi[,,t+1])%*%(Md[,t+1,i,j]*qb)+(Pi[,,t+1])%*%(M[,t+1,i]*qbd)
			if(mod>1){
				if(t>=mod) qbd = (Pid[,,2,j])%*%(M[,t+1,i]*qb)+(Pi[,,t+1])%*%(Md[,t+1,i,j]*qb)+(Pi[,,t+1])%*%(M[,t+1,i]*qbd)
				if(t<mod) qbd = (Pid[,,1,j])%*%(M[,t+1,i]*qb)+(Pi[,,t+1])%*%(Md[,t+1,i,j]*qb)+(Pi[,,t+1])%*%(M[,t+1,i]*qbd)
			}
			Qbd[,t,i,j] = qbd	
		}
	}	
}


# posterior probabilities
F1d = array(0,c(k,TT,ns,lth))
F1t = Q*Qb
for(j in 1:lth){
	F1dj = Qd[,,,j]*Qb+Q*Qbd[,,,j]
	for(i in 1:ns) F1d[,,i,j] = F1dj[,,i]/pv[i]-F1t[,,i]*pvd[i,j]/pv[i]^2
}
F2d = array(0,c(k,k,TT,ns,lth))
for(j in 1:lth){
	QbMdj = Qbd[,,,j]*M+Qb*Md[,,,j]
	for(i in 1:ns) for(t in 2:TT){
		if(mod==0) F2d[,,t,i,j] = Pid[,,t,j]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]+Pi[,,t]*(Qd[,t-1,i,j]%o%QbM[,t,i])/pv[i]+
                 	   Pi[,,t]*(Q[,t-1,i]%o%QbMdj[,t,i])/pv[i]-Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])*pvd[i,j]/pv[i]^2  
		if(mod==1) F2d[,,t,i,j] = Pid[,,j]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]+Pi[,,t]*(Qd[,t-1,i,j]%o%QbM[,t,i])/pv[i]+
                 	   Pi[,,t]*(Q[,t-1,i]%o%QbMdj[,t,i])/pv[i]-Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])*pvd[i,j]/pv[i]^2                
		if(mod>1){
			if(t <=mod) F2d[,,t,i,j] = Pid[,,1,j]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]+Pi[,,t]*(Qd[,t-1,i,j]%o%QbM[,t,i])/pv[i]+
                 	   Pi[,,t]*(Q[,t-1,i]%o%QbMdj[,t,i])/pv[i]-Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])*pvd[i,j]/pv[i]^2 
            if(t>mod) F2d[,,t,i,j] = Pid[,,2,j]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]+Pi[,,t]*(Qd[,t-1,i,j]%o%QbM[,t,i])/pv[i]+
                 	   Pi[,,t]*(Q[,t-1,i]%o%QbMdj[,t,i])/pv[i]-Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])*pvd[i,j]/pv[i]^2 
		}
	} 
}

out = list(lk=lk,sc=sc,F1=F1,F2=F2,F1d=F1d,F2d=F2d,M=M,Q=Q,pv=pv)

}

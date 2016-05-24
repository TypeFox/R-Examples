mica <-
function(W0=0,n.b=1000, maxit = 200, maxN=75, l.b=-10, u.b=10, N0=19, epsilon=10^(-4), hc=0, alpha=1,ind=500,nproc=10){
	
	#	added following lines to avoid package notes
	files=get("files");
	X=get("X");
	N.s=get("N.s");
	m=get("m");
	n=get("n");
	
	load(files[1])
	x=X
	T=nrow(x)
	V=ncol(x)
	abs.vec=c()

	## Find the binning parameters		
	counts=matrix(0,T,n.b)
	Tj=matrix(0,T,n.b)	
	
	for(subj in 1:N.s){
		load(files[subj])
		x=X
		tr=trans(t(W0[[subj]]%*%x),t(W0[[subj]]))
		W0[[subj]]=t(tr[[2]])
	}
	
	S0=t(tr[[1]])
	S.temp=pmax(pmin(S0,u.b),l.b)
	fout=fbins(t(S.temp),T,n.b)
	bins=fout[[1]]
	Tj=fout[[2]]	
	for(k in 1:T){
		counts[k,]=table(cut(S.temp[k,], breaks = bins[k,],include.lowest=TRUE))	}

	counts[counts==0]=1
	Co=counts/V

	### Create lists to include the means and weights of the mixture density
	mu.list=c()
	mu.full.list=c()
	theta.list=c()

	### Vector of starting values of variance components and mean vectors
	sigma.vec=rep(1,m)
	for(k in 1:m){
		sigma.vec[k]=2*(Tj[k,n.b-1]-Tj[k,2])/(3*(N0-1))
		mu.list=c(mu.list,list(seq(Tj[k,2],Tj[k,n.b-1],l=N0)))
		mu.full.list=c(mu.full.list,list(seq(Tj[k,2],Tj[k,n.b-1],l=2*N0-1)))
	}
	N=N0
	n.it=1
	abs.diff=1

	while((abs.diff>epsilon)&(N<=min(n-1,maxN))){
		j=1
		
		while((abs.diff>epsilon)&(2*j<(length(mu.full.list[[1]])/2))&(N<=min(V-1,maxN))){
	
			### Estimate the weights of the densities by using the EM algorithm
			theta.list=c()
			for(k in 1:m){	
				#### Estimate the density for this sample
				theta.hat=weight.seq(Tj[k,],Co[k,],mu.list[[k]],sigma.vec[k],0)
				theta.list=c(theta.list,list(theta.hat))
				}

			### update W0 using the computed mixture weights
			W1=c()
		 	inc=LW.grad(W0,theta.list,mu.list,sigma.vec,ind=ind,nproc=nproc,n,m,X,l.b,u.b,N.s,files)
			
			for(subj in 1:N.s){	
				W0.v=matrix(t(W0[[subj]]),m^2,1)
				#w1=matrix(W0.v+alpha*inc[[subj]][[1]],m,m)
				w1=matrix(W0.v+alpha*inc[[subj]],m,m)
				W1=c(W1,list(t(w1)))
			}
			
				###Set the value of the hidden sources using the last value of W

			for(subj in 1:N.s){
				load(files[subj])
				x=X
				tr=trans(t(W1[[subj]]%*%x),t(W1[[subj]]))
				W1[[subj]]=t(tr[[2]])
				}
			
			S0=t(tr[[1]])
			print(max((abs(cor(t(S0)))-diag(rep(1,m)))))
			S.temp=pmax(pmin(S0,u.b),l.b)
			fout=fbins(t(S.temp),m,n.b)
			bins=fout[[1]]
			Tj=fout[[2]]	
			for(k in 1:m){
				counts[k,]=table(cut(S.temp[k,], breaks = bins[k,],include.lowest=TRUE))			}

			counts[counts==0]=1
			Co=counts/V
			
			abs.diff=max(abs(W1[[1]]-W0[[1]]))
			for(subj in 1:N.s){
			abs.diff=max(abs.diff,max(abs(W1[[subj]]-W0[[subj]])))}
			abs.vec=c(abs.vec,abs.diff)
			print(abs.diff)
			W0=W1

			###Set the value of the hidden sources using the last value of W
			
			### Add two means
			N=N+2*(N<=maxN)
			if((abs.diff>epsilon)&(N<=maxN)){
				for(k in 1:m){
					mu.list[[k]]=sort(c(mu.list[[k]],mu.full.list[[k]][2*j],mu.full.list[[k]][length(mu.full.list[[k]])-2*j+1]))
}}

			j=j+1
			n.it=n.it+1					
		}
		
		mu.full.list=c()
		for(k in 1:m){
			mu.full.list=c(mu.full.list,list(seq(min(mu.list[[k]]),max(mu.list[[k]]),l=2*N-1)))}
	}

	while((abs.diff>epsilon)&(n.it<=maxit)){
	
			theta.list=c()
			### Estimate the weights of the densities by using the EM algorithm
			for(k in 1:m){	
				#### Estimate the density for this sample
				theta.hat=weight.seq(Tj[k,],Co[k,],mu.list[[k]],sigma.vec[k],0)
				theta.list=c(theta.list,list(theta.hat))
			}

			### update W0 using the computed mixture weights
			W1=c()
		 	inc=LW.grad(W0,theta.list,mu.list,sigma.vec,ind=ind,nproc=nproc,n,m,X,l.b,u.b,N.s,files)
			
			for(subj in 1:N.s){	
				W0.v=matrix(t(W0[[subj]]),m^2,1)
				#w1=matrix(W0.v+alpha*inc[[subj]][[1]],m,m)
				w1=matrix(W0.v+alpha*inc[[subj]],m,m)
				W1=c(W1,list(t(w1)))
			}
			
				###Set the value of the hidden sources using the last value of W

			for(subj in 1:N.s){
				load(files[subj])
				x=X
				tr=trans(t(W1[[subj]]%*%x),t(W1[[subj]]))
				W1[[subj]]=t(tr[[2]])
				}
			S0=t(tr[[1]])
			print(max((abs(cor(t(S0)))-diag(rep(1,m)))))
			S.temp=pmax(pmin(S0,u.b),l.b)
			fout=fbins(t(S.temp),m,n.b)
			bins=fout[[1]]
			Tj=fout[[2]]	
			for(k in 1:m){
				counts[k,]=table(cut(S.temp[k,], breaks = bins[k,],include.lowest=TRUE))			}

			counts[counts==0]=1
			Co=counts/V
							
			abs.diff=max(abs(W1[[1]]-W0[[1]]))
			for(subj in 1:N.s){
			abs.diff=max(abs.diff,max(abs(W1[[subj]]-W0[[subj]])))}
			abs.vec=c(abs.vec,abs.diff)
			print(abs.diff)
			W0=W1
		
			n.it=n.it+1
		}
		
	ret=c(list(W1),list(c(N,n.it)),list(theta.list),list(mu.list),list(sigma.vec),list(abs.vec))
	return(ret)
}

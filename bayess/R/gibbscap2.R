gibbscap2=function(nsimu,z){

# GIBBS SAMPLING FOR THE ARNASON-SCHWARZ CAPTURE-RECAPTURE MODEL

m=max(z)
T=dim(z)[2]
n=dim(z)[1]
p=array(0,c(nsimu,m))
phi=array(0,c(nsimu,m))
psi=array(0,c(m,m,nsimu))
latent=z

for (i in 1:n){

  for (t in 1:T){
	if (z[i,t]==0 & sum(z[i,t:T])!=0) 
		latent[i,t]=sample(1:m,1,prob=rep(1,m))
	if (z[i,t]==0 & sum(z[i,t:T])==0) 
		latent[i,t]=sample(1:(m+1),1,prob=c(rep(1,m),m))
	if (t!=1) if (latent[i,t-1]==m+1) latent[i,t]=m+1
	}
}


latentmean=latent
omega=rep(0,m+1)

for (s in 2:nsimu){

  for (r1 in 1:m) { for (r2 in 1:(m+1)){

         omega[r2]=sum(latent[,1:(T-1)]==r1 & latent[,2:T]==r2)
         }

  u=sum(z!=0 & latent==r1)
  v=sum(z==0 & latent==r1)
  p[s,r1]=rbeta(1,1+u,1+v)
  phi[s,r1]=rbeta(1,1+sum(omega[1:m]),1+omega[m+1])
  psi[r1,,s]=rdirichlet(1,rep(1,m)+omega[1:m])
  }

  tt=matrix(rep(phi[s,],m),m,byrow=T)
  q=rbind(tt*psi[,,s],rep(0,m))
  q=cbind(q,1-apply(q,1,sum))

  for (i in 1:n){

    if (z[i,1]==0) latent[i,1]=sample(1:(m+1),1,prob=q[,latent[i,2]]*(1-c(p[s,],0)))

    for (t in 2:(T-1)){

       if (z[i,t]==0) latent[i,t]=sample(1:(m+1),1,prob=q[latent[i,t-1],]*q[,latent[i,t+1]]*(1-c(p[s,],0)))
       }

    if (z[i,T]==0) latent[i,T]=sample(1:(m+1),1,prob=q[latent[i,T-1],]*(1-c(p[s,],0)))
    }

    latentmean=latentmean+latent
  }

latentmean=latentmean/nsimu
list(p=p,phi=phi,psi=psi,late=latentmean/nsimu)
}

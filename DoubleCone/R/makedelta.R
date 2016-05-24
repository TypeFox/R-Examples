makedelta <-
function(x,sh){
	n=length(x)
	xs=sort(x)
	xu=1:n*0
# find unique x values
	xs=round(xs*1e12)/1e12
	xu=unique(xs)
	n1=length(xu)
	sm=1e-8
# make bmat --equality constraints
	obs=1:n
	if(n1<n){
		bmat=matrix(0,nrow=n-n1,ncol=n)
		row=0
		for(i in 1:n1){
			cobs=obs[x==xu[i]]
			nr=length(cobs)
			if(nr>1){
				for(j in 2:nr){
					row=row+1
					bmat[row,cobs[1]]=-1;bmat[row,cobs[j]]=1
				}
			}
		}	
	}
# "hyper" convex or concave
	if(sh==3){
		amat=matrix(0,nrow=n1-3,ncol=n)
		for(i in 1:(n1-3)){
			c1=min(obs[abs(x-xu[i])<sm]);c2=min(obs[abs(x-xu[i+1])<sm])
			c3=min(obs[abs(x-xu[i+2])<sm]);c4=min(obs[abs(x-xu[i+3])<sm])
			amat[i,c1]=-(xu[i+3]-xu[i+2])*(xu[i+3]-xu[i+1])*(xu[i+2]-xu[i+1])
			amat[i,c2]=(xu[i+3]-xu[i])*(xu[i+3]-xu[i+2])*(xu[i+2]-xu[i])
			amat[i,c3]=-(xu[i+3]-xu[i])*(xu[i+3]-xu[i+1])*(xu[i+1]-xu[i])
			amat[i,c4]=(xu[i+1]-xu[i])*(xu[i+2]-xu[i])*(xu[i+2]-xu[i+1])
		}
	}else if(sh==2){
#  convex or concave
		amat=matrix(0,nrow=n1-2,ncol=n)
		for(i in 1:(n1-2)){
			c1=min(obs[abs(x-xu[i])<sm]);c2=min(obs[abs(x-xu[i+1])<sm]);c3=min(obs[abs(x-xu[i+2])<sm])
			amat[i,c1]=xu[i+2]-xu[i+1];amat[i,c2]=xu[i]-xu[i+2];amat[i,c3]=xu[i+1]-xu[i]
		}
#  increasing or decreasing
	}else{
		amat=matrix(0,nrow=n1-1,ncol=n)
		for(i in 1:(n1-1)){
			c1=min(obs[abs(x-xu[i])<sm]);c2=min(obs[abs(x-xu[i+1])<sm])
			amat[i,c1]=-1;amat[i,c2]=1
		}		
	}
	if(n1<n){
		wmat=matrix(0,nrow=n,ncol=n1)
		for(i in 1:n1){wmat[abs(x-xu[i])<sm,i]=1}
		atil=amat%*%wmat
		delta=t(wmat%*%t(atil)%*%solve(atil%*%t(atil)))
	}else{delta=solve(amat%*%t(amat))%*%amat}
	dr=length(delta)/n
	if(sh==3){
		pr2=cbind(1:n*0+1,x,x^2/max(x^2))
		prmat=pr2%*%solve(t(pr2)%*%pr2)%*%t(pr2)
		for(i in 1:dr){delta[i,]=delta[i,]-t(prmat%*%delta[i,])}
	}else if(sh==2){
		pr1=cbind(1:n*0+1,x)
		prmat=pr1%*%solve(t(pr1)%*%pr1)%*%t(pr1)
		for(i in 1:dr){delta[i,]=delta[i,]-t(prmat%*%delta[i,])}
	}else{
		for(i in 1:dr){delta[i,]=delta[i,]-mean(delta[i,])}
	}
	for(i in 1:dr){delta[i,]=delta[i,]/sqrt(sum(delta[i,]^2))}
	delta
}

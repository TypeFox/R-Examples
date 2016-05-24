fbins <-
function(s0,m,n.b){
	
	Tj=matrix(0,m,n.b)
	prob=seq(0,1,l=n.b)
	bins=matrix(0,m,n.b+1)
		
	for(k in 1:m){
		temp=quantile(s0[,k],prob=prob)
		temp1=temp
		ind=2
		while((ind<=n.b)&(temp1[ind]!=max(temp1))){
			ind=ind+1
			}
		if(ind<n.b){
		temp[(ind-1):n.b]=seq(temp1[ind-1],temp1[n.b],l=n.b-ind+2)}
		temp1=temp
		ind=2
		while(ind<=(n.b-1)){
			while((round(temp[ind],5)!=round(temp[ind-1],5))&(ind<=(n.b-1))){
				ind=ind+1
			}
			if((ind<n.b)&(round(temp[ind],5)==round(temp[ind-1],5))){
				ct=0
				while(((ind+ct)<n.b)&(round(temp[ind+ct],5)==round(temp[ind+ct-1],5))){ct=ct+1}
				
				temp1[(ind-1):(ind+ct+1)]=seq(temp[ind-1],temp[ind+ct+1],l=ct+3)
				ind=ind+ct
				}
			}
		Tj[k,]=temp1
	}
	
	for(k in 1:m){
		bins[k,1]=Tj[k,1]
		for(l in 2:n.b){
			delta1=Tj[k,l]-Tj[k,l-1]
			bins[k,l]=Tj[k,l]-1/2*delta1
		}
		bins[k,n.b+1]=Tj[k,n.b]
	}

	return(c(list(bins),list(Tj)))
	
	}

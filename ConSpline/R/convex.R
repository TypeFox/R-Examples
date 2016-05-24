convex <-
function(x,t){
	n=length(x)
	k=length(t)-2
	m=k+2
	sigma=matrix(1:m*n,nrow=m,ncol=n)
	dsigma=matrix(1:m*n,nrow=m,ncol=n)
	for(j in 1:(k-1)){
		i1=x<=t[j]
		sigma[j,i1] = 0
		dsigma[j,i1] = 0
	 	i2=x>t[j]&x<=t[j+1]
	 	sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
	 	dsigma[j,i2] = 3*(x[i2]-t[j])^2 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
	    i3=x>t[j+1]&x<=t[j+2]
	    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
	    dsigma[j,i3] = 1-3*(x[i3]-t[j+2])^2/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3
	    i4=x>t[j+2]
	    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
	    dsigma[j,i4]=1
	}
	i1=x<=t[k]
	sigma[k,i1] = 0
	dsigma[k,i1] = 0
	i2=x>t[k]&x<=t[k+1]
	sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
	dsigma[k,i2] = 3*(x[i2]-t[k])^2 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
	i3=x>t[k+1]
	sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
	dsigma[k,i3] = 1-3*(x[i3]-t[k+2])^2/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3
	i1=x<=t[2]
	sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3-(t[2]-t[1])^3/(t[2]-t[1])^2/3
	dsigma[k+1,i1]=1-3*(t[2]-x[i1])^2/(t[2]-t[1])^2/3
	i2=x>t[2]
	sigma[k+1,i2]=x[i2]-t[1]-(t[2]-t[1])^3/(t[2]-t[1])^2/3
	dsigma[k+1,i2]=1
	i1=x<=t[k+1]
	sigma[k+2,i1]=0
	dsigma[k+2,i1]=0
	i2=x>t[k+1]
	sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
	dsigma[k+2,i2]=3*(x[i2]-t[k+1])^2/(t[k+2]-t[k+1])^2/3
	
	for(i in 1:m){
		rng=max(sigma[i,])
		sigma[i,]=sigma[i,]/rng
		dsigma[i,]=dsigma[i,]/rng
	}
	
	ans=new.env()
	ans$sigma=sigma
	ans$dsigma=dsigma
	ans
}

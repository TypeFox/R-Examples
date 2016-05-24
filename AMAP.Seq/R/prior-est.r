
initial.uv=function(data){
	count=data$count
	size=data$size
	disp=data$dispersion
	treat=data$treat
	count[count==0]=min(count[count>0])/10
	m=count/size
	u=matrix(m[,treat==1],ncol=sum(treat==1))
	v=matrix(m[,treat==2],ncol=sum(treat==2))
    u=rowMeans(log(u))
    v=rowMeans(log(v))
	return(list(u=u,v=v))
}
prior.grid=function(u,v,n.grid=100){
	h=sqrt(sd(u)*sd(v))*1.06/length(u)^.2
	x=seq(min(u),max(u),length.out=n.grid)
	y=seq(min(v),max(v),length.out=n.grid)
    x=rep(x,n.grid)
    y=rep(y,each=n.grid)
    nE=length(u)
    w=rep(0,length(x))
    for(e in 1:nE)
    	w=w+dnorm((x-u[e])/h)*dnorm((y-v[e])/h)/h^2/nE
   # w[w==0]=min(w[w>0])*1e-10
	return(list(u=x,v=y,w=w))
}
posterior.uv=function(data,prior){
	count=data$count
	size=data$size
	disp=data$dispersion
	treat=data$treat
	nE=nrow(count)
    u=v=rep(0,nE)
	for(e in 1:nE){
		uv=posterior.uv.one(count[e,],size[e,],disp[e],treat,prior)
		u[e]=uv[1]
		v[e]=uv[2]
	}
	return(list(u=u,v=v))
}
posterior.uv.one=function(n,s,d,t,prior){
	u=prior$u
	v=prior$v
	w=prior$w
	f=log(w)
	for(i in 1:length(n)){
		if(t[i]==1){
			  m=s[i]*exp(u)
		}else m=s[i]*exp(v)	 
		f=f+lgamma(n[i]+1/d)-lgamma(n[i]+1)-lgamma(1/d)-log(1+m*d)/d-log(1+1/m/d)*n[i]
	}
	f=exp(f-max(f))
	u=sum(u*f)/sum(f)
	v=sum(v*f)/sum(f)
	return(c(u,v))	
}


prior.est=function(data,n.grid=100,n.iter=5){
	uv=initial.uv(data)
	prior0=list(u=uv$u,v=uv$v,w=rep(1/length(uv$u),length(uv$u)))
	uv.pdf.steps=list(prior0)
	for(i in 1:n.iter){
		prior=prior.grid(uv$u,uv$v,n.grid)
		uv=posterior.uv(data,prior)
		D=dist.CvM2D(prior0,prior)  #### Cramer von Mises distance
		#if(D<1e-2) break
		print(D)
		uv.pdf.steps[[i+1]]=prior
		prior0=prior
	}
	
	data$uv.pdf=prior
	data$u=uv$u
	data$v=uv$v
	data$uv.pdf.steps=uv.pdf.steps
	return(data)
}
dist.CvM2D=function(pdf1,pdf0){	
	cdf1=pdf2cdf(pdf1,pdf0)
    cdf0=pdf2cdf(pdf0,pdf0)
	d=sum((cdf0$f-cdf1$f)^2*pdf0$w)
	return(d)
}
pdf2cdf=function(pdf,grid=NULL){
	if(is.null(grid)) grid=pdf
	cdf=grid
	x=cdf[[1]]
	y=cdf[[2]]
	nZ=length(x)
	z=rep(0,nZ)
	for(i in 1:nZ){
		id=pdf[[1]]<=x[i] & pdf[[2]]<=y[i]
		z[i]=sum(pdf[[3]][id])
	}
	return(list(x=x,y=y,f=z/max(z)))
}











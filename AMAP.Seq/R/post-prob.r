

test.amap.as=function(data,prior,b1=NULL,b2=NULL,logfc=NULL){
	count=data$count
	size=data$size
	disp=data$dispersion
	treat=data$treat
	u=prior$u
	v=prior$v
	w=prior$w
	nE=nrow(count)
	if(is.null(b1)) b1=data$bg.u
	if(is.null(b2)) b2=data$bg.v
	if(is.null(logfc)) logfc=log(1.5)
	if(length(logfc)==1) logfc=rep(logfc,nE)
	P=B110=B111=B10=B01=B00=rep(0,nE)
	for(e in 1:nE){
		id=u>=b1[e] & v>=b2[e] & abs(u-v)<=logfc[e]
		B110[e]=int.nbinom.prior(count[e,],size[e,],disp[e],treat,u[id],v[id],w[id])
		
		id=u>=b1[e] & v>=b2[e] & abs(u-v)>logfc[e]	
		B111[e]=int.nbinom.prior(count[e,],size[e,],disp[e],treat,u[id],v[id],w[id])

		id=u>=b1[e] & v<b2[e]
		B10[e]=int.nbinom.prior(count[e,],size[e,],disp[e],treat,u[id],v[id],w[id])

		id=v>=b2[e] & u<b1[e]
		B01[e]=int.nbinom.prior(count[e,],size[e,],disp[e],treat,u[id],v[id],w[id])

		id=u<b1[e] & v<b2[e]
		B00[e]=int.nbinom.prior(count[e,],size[e,],disp[e],treat,u[id],v[id],w[id])
	}
	B11=LogSum(B110,B111)
        B=LogSum(B11,LogSum(B10,LogSum(B01,B00)))
	pr.fc=B110-B
	#pr.ex=LogSum(B00,LogSum(B10,B01))-B
	pr.ex=B00-B
        pr.as=LogSum(B00,B11)-B
	pr.as1=LogSum(B01,LogSum(B11,B00))-B
	pr.as2=LogSum(B10,LogSum(B11,B00))-B
	pvalues=data.frame(logFC=pr.fc,Expr=pr.ex,AS=pr.as,AS1=pr.as1,AS2=pr.as2)
	return(pvalues)
}
test.amap.gene=function(data,prior,FC=NULL){
	count=data$count
	size=data$size
	disp=data$dispersion
	treat=data$treat
	u=prior$u
	v=prior$v
	w=prior$w
	nE=nrow(count)
	if(is.null(FC)) FC=1
        if(FC==1) FC=1.001
        logfc=log(FC)
	if(length(logfc)==1) logfc=rep(logfc,nE)
	B1=B0=rep(0,nE)
	for(e in 1:nE){
		id= abs(u-v)<=logfc[e]
		B0[e]=int.nbinom.prior(count[e,],size[e,],disp[e],treat,u[id],v[id],w[id])
		
		id= abs(u-v)>logfc[e]	
		B1[e]=int.nbinom.prior(count[e,],size[e,],disp[e],treat,u[id],v[id],w[id])

	}
	B=LogSum(B0,B1)
	log.pr=B0-B
	pvalues=data.frame(log.post.prob=log.pr,post.prob=exp(log.pr))
	return(pvalues)
}

int.nbinom.prior=function(n,s,d,t,u,v,w){
	if(is.null(u)|is.null(v)) return(-Inf)
	f=log(w)
	for(i in 1:length(n)){
		if(t[i]==1){
			  m=s[i]*exp(u)
		}else m=s[i]*exp(v)	 
		f=f+lgamma(n[i]+1/d)-lgamma(n[i]+1)-lgamma(1/d)-log(1+m*d)/d-log(1+1/m/d)*n[i]
	}
	f=max(f)+LogSum(exp(f-max(f)))
	return(f)
}
LogSum=function(a,b){
	maxab=a
	maxab[a<b]=b[a<b]
	maxab[maxab==-Inf]=0
	ab=log(exp(a-maxab)+exp(b-maxab))+maxab
	return(ab)
}




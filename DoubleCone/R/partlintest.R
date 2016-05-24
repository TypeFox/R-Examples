partlintest <-
function(x,y,zmat,h0=0,nsim=1000){
	n=length(y)
	if(length(x)!=n){print("ERROR: x and y must be vectors of the same length")}
	if(length(zmat)==n){zmat=matrix(zmat,ncol=1)}
	if(!is.matrix(zmat)){print("ERROR: zmat must be an n by k matrix")}
	delta=makedelta(x,h0+1)
	m=dim(delta)[1]
	one=1:n*0+1
	if(max(abs(one-zmat%*%solve(t(zmat)%*%zmat)%*%t(zmat)%*%one))>1e-8){
		zmat=cbind(one,zmat)
	}
	k=dim(zmat)[2]
	if(h0==0){
		vmat=zmat
	}else if(h0==1){
		vmat=cbind(zmat,x)
	}else if(h0==2){
		vmat=cbind(zmat,x,x^2/max(x^2))
	}else{
		vmat=zmat
	}
	pmat=vmat%*%solve(t(vmat)%*%vmat)%*%t(vmat)
	p0k=pmat%*%y
	dtil=delta
	for(i in 1:m){dtil[i,]=delta[i,]-pmat%*%delta[i,]}
	ans1=coneB(y,dtil)
	p1k=ans1$yhat+p0k
	ans2=coneB(y,-dtil)
	p2k=ans2$yhat+p0k
	s0=sum((y-p0k)^2)
	s1=sum((y-p1k)^2)
	s2=sum((y-p2k)^2)
	tstat=(s0-min(s1,s2))/s0
	if(nsim>0){
		tdist=1:nsim
		for(isim in 1:nsim){
			ysim=rnorm(n)
			p0=pmat%*%ysim
			ans1=coneB(ysim,dtil)
			p1=ans1$yhat+p0
			ans2=coneB(ysim,-dtil)
			p2=ans2$yhat+p0
			s0=sum((ysim-p0)^2)
			s1=sum((ysim-p1)^2)
			s2=sum((ysim-p2)^2)
			tdist[isim]=(s0-min(s1,s2))/s0
		}
		pval=sum(tdist>tstat)/nsim
	}
	ans=new.env()
	ans$p0=p0k
	ans$p1=p1k
	ans$p2=p2k
	ans$pval=pval
	ans
}

library('MASS')

#logit
logit=function(a){
	return(qlogis(a))
}

#sign of estimator beta(exclude intercept)
sgn=function(X, y, y.mean,n,p,beta=0, gamma=0, lambda=0){
	sgn.cor=(1/n)*t(X[,-(p+1)])%*%(y-y.mean)
	if(gamma!=0){
		sgn.exc=(1:p)[abs(beta)>=lambda*gamma]
		sgn.cor[sgn.exc]=0
	}
	return(sign(sgn.cor))
		
}

#second order derivative of likelihood
lik.cov=function(X, y.mean, family, index,n,p){
	Xsub=X[,c(index,p+1)]
	y.mean=as.vector(y.mean)
	if(family=='poisson') cov.mat=diag(y.mean, nrow=length(y.mean))
	if(family=='binomial') cov.mat=diag(y.mean*(1-y.mean),nrow=length(y.mean))
	return((1/n)*t(Xsub)%*%cov.mat%*%Xsub)
}

#derivative of beta (exclude intercept) with respect to lambda
deri=function(cov,sgn,index,beta=0,gamma=0,lambda=0,least.eigen=0){
	sgn=sgn[index]
	dim.cov=dim(cov)[1]
	temp.cov=cov[-dim.cov, -dim.cov]
	if(gamma==0) inv=ginv(temp.cov)
	else{
		sub.beta=beta[index]
		exc=(abs(sub.beta)<=lambda*gamma+.Machine$double.eps)
		deri.pen=least.eigen/gamma*diag(exc,nrow=length(exc))
		inv=ginv(temp.cov-deri.pen)
		sgn=exc*sgn
	}	
	return(-inv%*%sgn)
}

#least eigenvalue used in mcp case
least.eigen=function(cov){
	ei=eigen(cov)$values
	return(min(Re(ei),1))
}

#covariance matrix used in quadratic approximation	
W.matrix=function(X,beta,family){
	eta=as.vector(exp(-X%*%beta))
	if(family=='poisson'){
		W=diag(1/eta)
	}
	if(family=='binomial'){
		W=diag(1/(1+eta)/(1+1/eta))
	}
	return(W)
}


#approximate response used in mcp case for correction
y.app=function(index, p, beta, X, y, family){
	index=c(index, p+1)
	Xsub=as.matrix(X[,index])
	sub.beta=beta[index]
	eta=as.vector(exp(-X%*%beta))
	if(family=='poisson'){
	    y.app.b=y-1/eta
	    W=diag(1/eta)	
	} 
	else{
        y.app.b=y-1/(1+eta)
        W=diag(1/(1+eta)/(1+1/eta))		
	} 
	
	return(Xsub%*%sub.beta+ginv(W)%*%y.app.b)
}


#first order derivative used in Newton Raphson update
NR.deri=function(X, y, lambda, n, index, y.mean){
	Xsub=X[,index]
	inter.deri=mean(y.mean-y)
	beta.deri=as.vector((1/n)*t(Xsub)%*%(y.mean-y))
	beta.shr.deri=sign(beta.deri)*max(abs(beta.deri)-lambda, 0)
	return(c(beta.shr.deri, inter.deri))
}

#second order derivative used in Newton Raphson update
NR.secderi=function(X, y.mean,n,p, index, family){
	Xsub=X[,c(index,p+1)]
	if(family=='binomial') var=(1/n)*t(Xsub)%*%diag(y.mean*(1-y.mean))%*%Xsub
	if(family=='poisson') var=(1/n)*t(Xsub)%*%diag(y.mean)%*%Xsub
	return(var)
}

#first order derivative used in Newton Raphson MCP case
mcp.NR.deri=function(X, y, n, index, beta,sgn, y.app,lambda,gamma, W, p){
	Xsub=as.matrix(X[,c(index, p+1)])
	sub.beta=beta[c(index, p+1)]
	mcp.deri.sgn=sgn[index]
	#mcp.deri.sgn=sign(sub.beta)
        eta=as.vector(exp(-X%*%beta))
	
	mcp.deri.f=as.vector((-1/n)*t(Xsub)%*%W%*%(y.app-Xsub%*%sub.beta))
	mcp.deri.d=pmax(lambda*(1-abs(sub.beta[-length(sub.beta)])/(lambda*gamma)), 0)*mcp.deri.sgn
	mcp.deri.d=c(mcp.deri.d, 0)
	return(mcp.deri.f+mcp.deri.d)
}

#first order derivative used in Newton Raphson MCP case
#mcp.NR.secderi=function(X, index, beta, lambda, gamma,n, W,least.eigen, p){
#	index=c(index, p+1)
#	Xsub=X[,index]
#	sub.beta=beta[index]
#	exc=(abs(sub.beta)<=lambda*gamma+.Machine$double.eps)
	
#	secderi=(1/n)*t(Xsub)%*%W%*%Xsub-least.eigen*diag(c(rep((1/gamma),length(index)-1),0)*exc, nrow=length(index))
#	return(secderi)
#}


mcp.NR.secderi=function(X, index, beta, lambda, gamma,n, W,least.eigen, p){
	index=c(index, p+1)
	n.ind=length(index)
	Xsub=X[,index]
	sub.beta=beta[index]
	exc=(abs(sub.beta)<=lambda*gamma+.Machine$double.eps)
	#cat('exc=  ',exc,'\n')
	pro.mat=(1/n)*t(Xsub)%*%W%*%Xsub
	pro.mat.e=eigen(pro.mat)
	vec=pro.mat.e$vectors
	val=pro.mat.e$values
	concave=rep(1, n.ind)-c(rep(1/gamma, n.ind-1),0)*exc
    ada=vec%*%diag(val*concave, nrow=n.ind)%*%solve(vec)
	secderi=ada
	return(secderi)
}


#weight used in quadratic approx. in CD algorithm
qua.weight=function(y.mean,family){
	if(family=='binomial') return(y.mean*(1-y.mean))
	if(family=='poisson') return(y.mean)
}

#response used in quadratic approx. in CD algorithm
qua.res=function(X,beta,y,y.mean,family){
	if(family=='binomial') return(X%*%beta+(y-y.mean)/(y.mean*(1-y.mean)))
	if(family=='poisson') return(X%*%beta+(y-y.mean)/(y.mean))
}


#soft threshold
soft.thre=function(z, lambda){
	return(sign(z)*max(abs(z)-lambda, 0))
}

#negative log-liklihood
loglik=function(X, y, beta, family){
	link=as.vector(X%*%beta)
	if(family=='poisson') return(sum(exp(link)-y*link))
	if(family=='binomial') return(sum(log(1+exp(link))-y*link))
}

#bic
bic=function(n, loglik, beta){
	return(2*loglik+sum(beta!=0)*log(n))
}

#aic
aic=function(loglik, beta){
	return(loglik+sum(beta!=0))
}

#ebic
ebic=function(loglik, beta, n, p){
	K=sum(beta!=0)
	return(2*loglik+K*log(n)+2*log(choose(p,K)))
}


sub.vec=function(a,b){
	if(length(b)>0){
		for(i in 1:length(b)){
			if(length((1:length(a))[(a==b[i])])>0) a=a[-(1:length(a))[(a==b[i])]]
		}
	}
	if(length(a)>0) return(a)
	else return(TRUE)
}

#cross validation criterion
loss = function(y.hat,y.sub,family){
	if(family=='poisson'){
		deveta=y.sub*log(y.hat)-y.hat
		devy=y.sub*log(y.sub)-y.sub
		devy[y.sub==0]=0
		return(2*apply((devy-deveta),2,mean))
		}
		
	if(family=='binomial'){
		if(sum(y.sub==1)==1) part.1=sum(log(y.hat[y.sub==1,]))
	    else part.1=apply(log(y.hat[y.sub==1,]),2,sum)
	    if(sum(y.sub==0)==1) part.2=sum(log(1-y.hat[y.sub==0,]))
	    else part.2=apply(log(1-y.hat[y.sub==0,]),2,sum)
		return((part.1+part.2)/(-1*length(y.sub)))
		}
	}


#CDL = function(X, y, beta, family, max.iter, lambda, eps ){
#	dyn.load('CDL.so')
#	
#	n=nrow(X)
#	p=ncol(X)
#	
#	index = (1:(p-1))[beta[1:(p-1)]!=0]
#	p0=length(index)
#	if(family=='binomial') fam=2
#	if(family=='poisson') fam=1
#	paraIn = c(eps, max.iter, lambda)
#        bI = 0
#	
#	out = .Fortran('CDL',
#	       X = as.double(X), y=as.double(y), beta=as.double(beta),
#	       n = as.integer(n), p=as.integer(p), bI=as.integer(bI),
#	       fam=as.integer(fam), maxIte=as.integer(max.iter), p0=as.integer(p0),
#	       ind=as.integer(index), lam=as.double(lambda), eps=as.double(eps), paraIn=as.double(paraIn)
#	)
#
#  	a=list(beta=out$beta, break.ind=out$bI)
#    return(a)
#}

#CDM = function(X, y, beta, family, max.iter, index, lambda, eps, gamma){
#	dyn.load('CDM.so')
#	
#	n=nrow(X)
#	p=ncol(X)
#	#index = (1:(p-1))[beta[1:(p-1)]!=0]
#        p0=length(index)
#	
#	if(family=='binomial') fam=2
#	if(family=='poisson') fam=1
#	paraIn = c(eps, max.iter, lambda, gamma)
#	
#	out = .Fortran('CDM',
#	       X=as.double(X), y=as.double(y), p=as.integer(p), 
#               n=as.integer(n), beta=as.double(beta),
#               paraIn=as.double(paraIn), fam=as.integer(fam), 
#	       bI=as.integer(0),
#	       ind=as.integer(index), p0=as.integer(p0)
#	)
#	a=list(beta=out$beta, break.ind=out$bI) 
#	return(a)
#}

cd_lasso1 = function(X, y, beta, epsilon, max.iter, lambda, family, bInd){
	#dyn.load('cd_lasso1.so')
	
	p = ncol(X)
	n = nrow(X)
	para.in = c(epsilon, max.iter)
	bInd = bInd
	if(family=='binomial') family=2
	if(family=='poisson') family=1
	out = .Fortran('cd_lasso1',
	       X = as.double(X),
	       y = as.double(y),
	       p = as.integer(p),
	       n = as.integer(n),
	       beta = as.double(beta),
	       paraIn = as.double(para.in),
	       family = as.integer(family),
	       bInd = as.logical(bInd),
	       lam = as.double(lambda)
	       )
	       	       
return(list(beta=out$beta, bInd=out$bInd))		
}

#cd_lasso2 = function(X, y, beta, epsilon, max.iter, lambda, family, bInd){
#	dyn.load('cd_lasso2.so')
#	
#	p = ncol(X)
#	n = nrow(X)
#	para.in = c(epsilon, max.iter)
#	bInd = bInd
#	if(family=='binomial') family=2
#	if(family=='poisson') family=1
#	out = .Fortran('cd_lasso2',
#	       X = as.double(X),
#	       y = as.double(y),
#	       p = as.integer(p),
#	       n = as.integer(n),
#	       beta = as.double(beta),
#	       paraIn = as.double(para.in),
#	       family = as.integer(family),
#	       bInd = as.logical(bInd),
#	       lam = as.double(lambda)
#	       )
#	       	       
#return(list(beta=out$beta, bInd=out$bInd))		
#}

cd_mcp1 = function(X, y, beta, epsilon, max.iter, lambda, gamma, family, bInd){
	#dyn.load('cd_mcp1.so')
	p = ncol(X)
	n = nrow(X)
	para.in = c(epsilon, max.iter, gamma)
	#print(para.in)
	bInd = bInd
	if(family=='binomial') family=2
	if(family=='poisson') family=1
	out = .Fortran('cd_mcp1',
	       X = as.double(X),
	       y = as.double(y),
	       p = as.integer(p),
	       n = as.integer(n),
	       beta = as.double(beta),
	       paraIn = as.double(para.in),
	       family = as.integer(family),
	       bInd = as.logical(bInd),
	       lam = as.double(lambda)
	       )
	       
return(list(beta=out$beta, bInd=out$bInd))		
}

#cd_mcp2 = function(X, y, beta, epsilon, max.iter, lambda, gamma, family, bInd){
#	dyn.load('cd_mcp2.so')
	
#	p = ncol(X)
#	n = nrow(X)
#	para.in = c(epsilon, max.iter, gamma)
#	#print(para.in)
#	bInd = bInd
#	if(family=='binomial') family=2
#	if(family=='poisson') family=1
#	out = .Fortran('cd_mcp2',
#	       X = as.double(X),
#	       y = as.double(y),
#	       p = as.integer(p),
#	       n = as.integer(n),
#	       beta = as.double(beta),
#	       paraIn = as.double(para.in),
#	       family = as.character(family),
#	       bInd = as.logical(bInd),
#	       lam = as.double(lambda)
#	       )
#	       
#return(list(beta=out$beta, bInd=out$bInd))		
#}

f.intercept = function(y, family){
	if(family == 'poisson') intercept = log(mean(y))
	if(family == 'binomial') intercept = log(mean(y)/(1-mean(y)))
	return(intercept)
}

#Lambdas are taken on log scale
nnlasso<-function(x,y,family=c("normal","binomial","poisson"),lambda=NULL,intercept=TRUE,normalize=TRUE,tau=1,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4,eps=1e-6,path=TRUE,SE=FALSE)
{
	if(missing(family)) family="normal"
  	switch(family, normal = nnlasso.normal(x,y,lambda,intercept,normalize,tau,tol,maxiter,nstep,min.lambda,eps,path,SE),binomial = nnlasso.binomial(x,y,lambda,intercept,normalize,tau,tol,maxiter,nstep,min.lambda,eps,path,SE), poisson = nnlasso.poisson(x,y,lambda,intercept,normalize,tau,tol,maxiter,nstep,min.lambda,eps,path,SE))
}
##############################################################################################################################################
nnlasso.normal<-function(x,y,lambda=NULL,intercept=TRUE,normalize=TRUE,tau=1,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4,eps=1e-6,path=TRUE,SE=FALSE)
{	
	np=dim(x)
	n=np[1]
	p=np[2]
	if (intercept)
	{
        	meanx = colMeans(x)
        	x = scale(x, meanx, FALSE)
        	meany = sum(y)/n
        	y = y - meany
    	} else {
        	meanx = rep(0, p)
        	meany = 0        	
    		}
    	if (normalize) 
  	{
        	normx = sqrt(colSums(x^2))
		x = scale(x, FALSE, normx)		
    	} else normx = rep(1, p)
	tx=t(x)
	xpy=tx%*%y
	max.lambda=max(abs(xpy))
	if (path)
	{
		stepsize=(log(min.lambda)-log(max.lambda))/(nstep-1)
		lambdas=exp(log(max.lambda)+stepsize*c(0:(nstep-1)))
	} else {
			nstep=2
			lambdas=c(max.lambda,lambda)
		}
	coef=matrix(0,nstep,p)	
	of.value=rep(0,nstep)
	coef[1,]=1e-2
	xbeta.old=x%*%coef[1,]
	if (n>p)
	{
		xpx=tx%*%x
		xpxbetaold=xpx%*%coef[1,]
	} else {
			xpx=NULL
			xpxbetaold=tx%*%xbeta.old
		}
	of.value[1]=-sum((y-xbeta.old)^2)/2-max.lambda*(tau*sum(coef[1,])+(1-tau)*sum(coef[1,]^2))
	lambda.iter=rep(0,nstep)
	g1=xpy
	for(iter in 2:nstep)
	{
		kkt=FALSE
		while(kkt==FALSE)
		{
			res=nnlasso.normal.lambda(n,p,x,y,xpx,xpy,beta.old=coef[iter-1,],tau,lambda1=lambdas[iter],tol,maxiter,xbeta.old,eps,SE)
			if (res$conv=="yes")
			{	
				coef[iter,]=res$beta.new
				xbeta=res$xbeta.new
				if (n>p)
				{
					xpxbetaold=xpx%*%res$beta.new
				} else xpxbetaold=tx%*%xbeta
				g1=xpy-xpxbetaold-lambdas[iter]*2*(1-tau)*coef[iter,]
				indices=NULL
				if (length(indices)==0)
				{
					xbeta.old=res$xbeta.new
					lambda.iter[iter]=res$iter
					of.value[iter]=res$ofv.new
					kkt=TRUE
				}
			} else stop("The algorithm did not converge")
		}
	}
	coef=scale(coef,center=FALSE,scale=normx)
	if(SE)
	{
		vcov=res$vcov
		vcov=vcov/normx	
		vcov=t(vcov)/normx
		se=sqrt(diag(vcov))
		if(intercept) 
		{
			se0=sqrt(t(meanx)%*%vcov%*%meanx)
			se=c(se0,se)
		}				
	} else se=NULL
	if (intercept) beta0=rep(meany,nstep)-coef%*%meanx else beta0=rep(0,nstep)
	L1norm=rowSums(abs(coef))
	norm.frac=L1norm/max(L1norm)
	obj=list(beta0=beta0,coef=coef,lambdas=lambdas,L1norm=L1norm,norm.frac=norm.frac,lambda.iter=lambda.iter,of.value=of.value,normx=normx,se=se)
	class(obj)='nnlasso'
	return(obj)
}
##############################################################################################################################################################
nnlasso.binomial<-function(x,y,lambda=NULL,intercept=TRUE,normalize=TRUE,tau=1,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4,eps=1e-6,path=TRUE,SE=FALSE)
{
	np=dim(x)
	n=np[1]
	p=np[2]	   	
    	if (intercept) 
	{
        	meanx = colMeans(x)
        	x = scale(x, meanx, FALSE)        	
    	} else {
        	meanx = rep(0, p)
        	meany = 0        	
    		}
    	if (normalize) 
  	{
        	normx = sqrt(colSums(x^2))           
        	x = scale(x, FALSE, normx)
    	} else normx = rep(1, p)
	tx=t(x)
	tx2=tx^2
	sumy=sum(y)
	g1=tx%*%(y-0.5)
	max.lambda=max(abs(g1))
	if (path)
	{
		stepsize=(log(min.lambda)-log(max.lambda))/(nstep-1)
		lambdas=exp(log(max.lambda)+stepsize*c(0:(nstep-1)))
	} else {
			nstep=2
			lambdas=c(max.lambda,lambda)
		}
	beta0=rep(0,nstep)
	coef=matrix(0,nstep,p)
	coef[1,]=1e-2
	xbeta.old=matrix(0,nstep,n)
	xbeta.old[1,]=x%*%coef[1,]
	Mu=matrix(0.5,nstep,n)
	lambda.iter=rep(0,nstep)
	of.value=rep(0,nstep)
	expeta.old=exp(xbeta.old[1,])
	of.value[1]=sum(y*xbeta.old[1,])-sum(log(1+expeta.old))-max.lambda*(tau*sum(coef[1,])+(1-tau)*sum(coef[1,]^2))
	Mu[1,]=expeta.old/(1+expeta.old)
	for(iter in 2:nstep)
	{
		dxkx0=ifelse(intercept,sum(Mu[iter-1,]*(1-Mu[iter-1,])),Inf)
		kkt=FALSE
		while(kkt==FALSE)
		{
			res=nnlasso.binomial.lambda(n,p,sumy,beta0.old=beta0[iter-1],beta1.old=coef[iter-1,],x,y,dxkx0,tau,lambda1=lambdas[iter],tol,maxiter,xbeta.old[iter-1,],mu1=Mu[iter-1,],eps,SE)
			if (res$conv=="yes") 
			{
				coef[iter,]=res$beta1.new
				g1=tx%*%(y-res$mu1)-2*lambdas[iter]*(1-tau)*coef[iter-1,]	
				indices=NULL
				if (length(indices)==0)
				{
					lambda.iter[iter]=res$iter
					of.value[iter]=res$ofv.new
					xbeta.old[iter,]=res$xbeta.new
					beta0[iter]=res$beta0.new
					Mu[iter,]=res$mu1
					kkt=TRUE
				}
			} else stop("The algorithm did not converge")	
		}
	}	
	coef=scale(coef,center=F,scale=normx)
	if(SE)
	{
		if(intercept)
		{
			vcov=res$vcov[2:(p+1),2:(p+1)]
			vcov=vcov/normx	
			vcov=t(vcov)/normx
			se=sqrt(diag(vcov))
			se0=sqrt(res$vcov[1,1]-t(meanx)%*%vcov%*%meanx)
			se=c(se0,se)
		} else {
				vcov=res$vcov
				vcov=vcov/normx	
				vcov=t(vcov)/normx
				se=sqrt(diag(vcov))
			}		
						
	} else se=NULL
	L1norm=rowSums(abs(coef))
	norm.frac=L1norm/max(L1norm)
	if (intercept) beta0=beta0-coef%*%meanx
	obj=list(beta0=beta0,coef=coef,lambdas=lambdas,L1norm=L1norm,norm.frac=norm.frac,lambda.iter=lambda.iter,of.value=of.value,normx=normx,se=se)
	class(obj)='nnlasso'
	return(obj)
}
##################################################################################################################################################################
nnlasso.poisson<-function(x,y,lambda=NULL,intercept=TRUE,normalize=TRUE,tau=1,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4,eps=1e-6,path=TRUE,SE=FALSE)
{
	np=dim(x)
	n=np[1]
	p=np[2]	   	
    	if (intercept) 
	{
        	meanx = colMeans(x)
        	x = scale(x, meanx, FALSE)        	
    	} else {
       		meanx = rep(0, p)
        	meany = 0        	
    		}
    	if (normalize) 
  	{
        	normx = sqrt(colSums(x^2))           
        	x = scale(x, FALSE, normx)
    	} else normx = rep(1, p)
	tx=t(x)
	tx2=tx^2
	sumy=sum(y)
	g1=tx%*%(y-1)
	max.lambda=max(abs(g1))
	if (path)
	{
		stepsize=(log(min.lambda)-log(max.lambda))/(nstep-1)
		lambdas=exp(log(max.lambda)+stepsize*c(0:(nstep-1)))
	} else {
			nstep=2
			lambdas=c(max.lambda,lambda)
		}
	beta0=rep(0,nstep)
	coef=matrix(0,nstep,p)
	coef[1,]=1e-2
	xbeta.old=matrix(0,nstep,n)
	xbeta.old[1,]=x%*%coef[1,]	
	lambda.iter=rep(0,nstep)
	Mu=matrix(1,nstep,n)
	Mu[1,]=exp(xbeta.old[1,])
	of.value=rep(0,nstep)
	of.value[1]=sum(y*xbeta.old[1,])-sum(Mu[1,])-max.lambda*(tau*sum(coef[1,])+(1-tau)*sum(coef[1,]^2))	
	for(iter in 2:nstep)
	{
		dxkx0=ifelse(intercept,sum(Mu[iter-1,]),Inf)
		kkt=FALSE
		while(kkt==FALSE)
		{
			res=nnlasso.poisson.lambda(n,p,sumy,beta0.old=beta0[iter-1],beta1.old=coef[iter-1,],x,y,dxkx0,tau,lambda1=lambdas[iter],tol,maxiter,xbeta.old[iter-1,],mu1=Mu[iter-1,],eps,SE)
			if (res$conv=="yes") 
			{
				coef[iter,]=res$beta1.new
				g1=tx%*%(y-res$mu1)-2*lambdas[iter]*(1-tau)*coef[iter-1,]	
				indices=NULL
				if (length(indices)==0)
				{
					lambda.iter[iter]=res$iter
					of.value[iter]=res$ofv.new
					beta0[iter]=res$beta0.new
					xbeta.old[iter,]=res$xbeta.new
					Mu[iter,]=res$mu1
					kkt=TRUE
				} 
			} else stop("The algorithm did not converge")	
		}
	}
	coef=scale(coef,center=F,scale=normx)
	if(SE)
	{
		if(intercept)
		{
			vcov=res$vcov[2:(p+1),2:(p+1)]
			vcov=vcov/normx	
			vcov=t(vcov)/normx
			se=sqrt(diag(vcov))
			se0=sqrt(res$vcov[1,1]-t(meanx)%*%vcov%*%meanx)
			se=c(se0,se)
		} else {
				vcov=res$vcov
				vcov=vcov/normx	
				vcov=t(vcov)/normx
				se=sqrt(diag(vcov))
			}		
						
	} else se=NULL
	L1norm=rowSums(abs(coef))
	norm.frac=L1norm/max(L1norm)
	if (intercept) beta0=beta0-coef%*%meanx
	obj=list(beta0=beta0,coef=coef,lambdas=lambdas,L1norm=L1norm,norm.frac=norm.frac,lambda.iter=lambda.iter,of.value=of.value,normx=normx,se=se)
	class(obj)='nnlasso'
	return(obj)
}
#################################################################################################################################################################
nnlasso.normal.lambda<-function(n,p,x,y,xpx,xpy,beta.old,tau,lambda1,tol,maxiter,xbeta.old,eps,SE=FALSE)
{
	epp=0.001
	if (n<=p) tx=t(x)
	ofv.old=-sum((y-xbeta.old)^2)/2-lambda1*(tau*sum(beta.old)+(1-tau)*sum(beta.old^2))
	for (iter in 1:maxiter)
	{
		if (n>p)
		{
			xpxbetaold=xpx%*%beta.old
		} else xpxbetaold=tx%*%xbeta.old
		g1minus<-g1<-xpy-xpxbetaold
		g1minus[which(g1>0)]=0
		b=beta.old*(g1-lambda1*tau-2*lambda1*(1-tau)*beta.old)/(lambda1*tau+2*lambda1*(1-tau)*beta.old-g1minus+1e-16)
		beta.new=beta.old+b
		xb=x%*%b
		xbeta.new=xbeta.old+xb					       
		ofv.new=-sum((y-xbeta.new)^2)/2-lambda1*(tau*sum(beta.new)+(1-tau)*sum(beta.new^2))
		delta=1
		t1=epp*(sum(g1*b))		
		while (ofv.new-delta*t1<ofv.old & delta>1e-5)			
		{	
			delta=delta/2
			beta.new=beta.old+delta*b			
			xbeta.new=xbeta.old+delta*xb
			ofv.new=-sum((y-xbeta.new)^2)/2-lambda1*(tau*sum(beta.new)+(1-tau)*sum(beta.new^2))
		}
		if (ofv.new-delta*t1<ofv.old & delta<=1e-5) 
		{
			beta.new=beta.old
			ofv.new=ofv.old	
			xbeta.new=xbeta.old
			break
		}
		if(abs(ofv.old-ofv.new)<=tol) break
		beta.old=beta.new		
		xbeta.old=xbeta.new
		ofv.old=ofv.new
	}
	if (iter<maxiter) conv="yes" else conv="no"
	if (SE==TRUE & conv=="yes")
	{
		grad=xpy-t(x)%*%xbeta.new-lambda1*tau-2*lambda1*(1-tau)*beta.new
    		Finv<-F<-xpx+2*lambda1*(1-tau) 
		index=which(beta.new<=eps & grad< -1e-2)
		if (length(index)>0) 
		{
			temp=F[-index,-index]
			tempinv=solve(temp)
			indexc=setdiff(1:p,index)
			if (length(indexc)>0) Finv[indexc,indexc]=tempinv
			Finv[index,index]=0
		} else Finv=solve(F)
		vc=Finv%*%F
		vcov=vc%*%t(Finv)			
	} else vcov=NULL
	res=list(beta.new=beta.new,conv=conv,iter=iter,ofv.new=ofv.new,xbeta.new=xbeta.new,vcov=vcov)
	return(res)		
}
####################################################################################################################################################
nnlasso.binomial.lambda<-function(n,p,sumy,beta0.old,beta1.old,x,y,dxkx0,tau,lambda1,tol,maxiter,xbeta.old,mu1,eps,SE=FALSE)
{
	epp=0.001
	tx=t(x)
	expeta.old=exp(xbeta.old)
	ofv.old=sum(y*xbeta.old)-sum(log(1+expeta.old))-lambda1*(tau*sum(beta1.old)+(1-tau)*sum(beta1.old^2))
	for (iter in 1:maxiter)
	{
		g0=sumy-sum(mu1)
		b0=g0/dxkx0
		beta0.new=beta0.old+b0
		ldashminus<-ldash<-tx%*%(y-mu1)
		ldashminus[which(ldash>0)]=0
		g1=ldash-lambda1*tau-2*lambda1*(1-tau)*beta1.old
		b1=beta1.old*g1/(lambda1*tau+2*lambda1*(1-tau)*beta1.old-ldashminus)
		beta1.new=beta1.old+b1
		t1=epp*(sum(g1*b1)+g0*b0)
		xb=x%*%b1+b0
		xbeta.new=xbeta.old+xb
		expeta.new=exp(xbeta.new)
		ofv.new=sum(y*xbeta.new)-sum(log(1+expeta.new))-lambda1*(tau*sum(beta1.new)+(1-tau)*sum(beta1.new^2))
		delta=1
		while (ofv.new-delta*t1<ofv.old & delta>1e-5)			
		{	
			delta=delta/2
			beta0.new=beta0.old+delta*b0 		
			beta1.new=beta1.old+delta*b1
			xbeta.new=xbeta.old+delta*xb
			expeta.new=exp(xbeta.new)
			ofv.new=sum(y*xbeta.new)-sum(log(1+expeta.new))-lambda1*(tau*sum(beta1.new)+(1-tau)*sum(beta1.new^2))
		}
		if (ofv.new-delta*t1<ofv.old & delta<1e-5) 
		{
			beta0.new=beta0.old										
			beta1.new=beta1.old										
			xbeta.new=xbeta.old
			ofv.new=ofv.old						
			break	
		}
		mu1=expeta.new/(1+expeta.new)									
		if(abs(ofv.old-ofv.new)<tol) break
		beta0.old=beta0.new
		beta1.old=beta1.new
		xbeta.old=xbeta.new
		ofv.old=ofv.new		
	}
	if (iter<maxiter) conv="yes" else conv="no"	 
	if(SE==TRUE & conv=="yes")
	{
		if (dxkx0!=Inf) 
		{
			x=cbind(rep(1,n),x)  
			p=p+1
			grad=t(x)%*%(y-mu1)-lambda1*tau-2*lambda1*(1-tau)*c(0,beta1.new)
		} else grad=t(x)%*%(y-mu1)-lambda1*tau-2*lambda1*(1-tau)*beta1.new		
    		k=diag(as.vector(mu1*(1-mu1)))
		xkx=t(x)%*%k%*%x
		Hinv<-H<-xkx+2*lambda1*(1-tau)
		if (dxkx0!=Inf) index=which(beta1.new<=eps & grad[2:p] < -1e-2) else index=which(beta1.new<=eps & grad < -1e-2)
		if (length(index)>0) 
		{
			if (dxkx0!=Inf) index=index+1
			temp=H[-index,-index]
			tempinv=solve(temp)
			indexc=setdiff(1:p,index)
			if (length(indexc)>0) Hinv[indexc,indexc]=tempinv
			Hinv[index,index]=0
		} else Hinv=solve(H)
		vc=Hinv%*%H
		vcov=vc%*%t(Hinv)		
	} else vcov=NULL
	res=list(beta0.new=beta0.new,beta1.new=beta1.new,conv=conv,iter=iter,ofv.new=ofv.new,xbeta.new=xbeta.new,mu1=mu1,vcov=vcov)
	return(res)
}
##############################################################################################################################################################
nnlasso.poisson.lambda<-function(n,p,sumy,beta0.old,beta1.old,x,y,dxkx0,tau,lambda1,tol,maxiter,xbeta.old,mu1,eps,SE=FALSE)
{
	epp=0.001
	tx=t(x)	
	ofv.old=sum(y*xbeta.old)-sum(mu1)-lambda1*(tau*sum(beta1.old)+(1-tau)*sum(beta1.old^2))
	for (iter in 1:maxiter)
	{
		g0=sumy-sum(mu1)
		b0=g0/dxkx0
		beta0.new=beta0.old+b0
		ldashminus<-ldash<-tx%*%(y-mu1)
		ldashminus[which(ldash>0)]=0
		g1=ldash-lambda1*tau-2*lambda1*(1-tau)*beta1.old
		b1=beta1.old*g1/(lambda1*tau+2*lambda1*(1-tau)*beta1.old-ldashminus)		
		beta1.new=beta1.old+b1
		t1=epp*(sum(g1*b1)+g0*b0)
		xb=x%*%b1+b0
		xbeta.new=xbeta.old+xb
		expeta.new=exp(xbeta.new)
		ofv.new=sum(y*xbeta.new)-sum(expeta.new)-lambda1*(tau*sum(beta1.new)+(1-tau)*sum(beta1.new^2))
		delta=1		
		while (ofv.new-delta*t1<ofv.old & delta>1e-5)			
		{	
			delta=delta/2	
			beta0.new=beta0.old+delta*b0		
			beta1.new=beta1.old+delta*b1
			xbeta.new=xbeta.old+delta*xb
			expeta.new=exp(xbeta.new)
			ofv.new=sum(y*xbeta.new)-sum(expeta.new)-lambda1*(tau*sum(beta1.new)+(1-tau)*sum(beta1.new^2))
		}
		if (ofv.new-delta*t1<ofv.old & delta<1e-5) 
		{
			beta0.new=beta0.old
			beta1.new=beta1.old										
			xbeta.new=xbeta.old
			ofv.new=ofv.old						
			break	
		}
		mu1=expeta.new									
		if(abs(ofv.old-ofv.new)<tol) break
		beta0.old=beta0.new
		beta1.old=beta1.new
		xbeta.old=xbeta.new
		ofv.old=ofv.new		
	}
	if (iter<maxiter) conv="yes" else conv="no"
	if(SE==TRUE & conv=="yes")
	{
		if (dxkx0!=Inf) 
		{
			x=cbind(rep(1,n),x)  
			p=p+1
			grad=t(x)%*%(y-mu1)-lambda1*tau-2*lambda1*(1-tau)*c(0,beta1.new)
		} else grad=t(x)%*%(y-mu1)-lambda1*tau-2*lambda1*(1-tau)*beta1.new
    		k=diag(as.vector(mu1))
		Hinv<-H<-t(x)%*%k%*%x+2*lambda1*(1-tau)
		if (dxkx0!=Inf) index=which(beta1.new<=eps & grad[2:p] < -1e-3) else index=which(beta1.new<=eps & grad < -1e-3)
		if (length(index)>0) 
		{
			if (dxkx0!=Inf) index=index+1
			temp=H[-index,-index]
			tempinv=solve(temp)
			indexc=setdiff(1:p,index)
			if (length(indexc)>0) Hinv[indexc,indexc]=tempinv
			Hinv[index,index]=0
		} else Hinv=solve(H)
		vc=Hinv%*%H
		vcov=vc%*%t(Hinv)	
	} else vcov=NULL
	res=list(beta0.new=beta0.new,beta1.new=beta1.new,conv=conv,iter=iter,ofv.new=ofv.new,xbeta.new=xbeta.new,mu1=mu1,vcov=vcov)
	return(res)
}
###############################################################################################################################################################################
predict.nnlasso<-function(object,mode=c("fraction","norm","lambda"),at=0,...)
{
	if (missing(mode)) mode="lambda"
	coef=object$coef
	beta0=object$beta0
	L1norm=object$L1norm
	if (mode=="lambda") values=object$lambdas else if (mode=="fraction") values=object$norm.frac else values=L1norm
	if (any(at>=max(values))) 
	{
		max.value=which(values==max(values))
		pred0=beta0[max.value]
		pred1=coef[max.value,] 
		
	} else if (any(at<=min(values))) 
	{
		min.value=which(values==min(values))
		pred0=beta0[min.value]
		pred1=coef[min.value,] 
		
	} else {
		value.up=min(values[values>at])
		value.down=max(values[values<at])
		uprownum=which(values==value.up)
		coef.up=coef[uprownum,]
		beta0.up=beta0[uprownum]
		downrownum=which(values==value.down)
		coef.down=coef[downrownum,]
		beta0.down=beta0[downrownum]
		if (value.up!=value.down) 
		{
			pred0=beta0.down+(beta0.up-beta0.down)*(at-value.down)/(value.up-value.down)
			pred1=coef.down +(coef.up-coef.down)*(at-value.down)/(value.up-value.down) 
		} else {
			pred0=beta0.up
			pred1=coef.up
			}
	      }
	pred=c(pred0,pred1)	
	return(pred)
}
########################################################################################################################################################################
plot.nnlasso<-function(x,xvar=c("lambda","L1norm","fraction of norm"),...)
{
	coef=x$coef
	lambda=x$lambdas
	norm.frac=x$norm.frac
	L1norm=x$L1norm	
	if (missing(xvar)) xvar="L1norm"
	if (xvar=="lambda") 
	{
		matplot(lambda,coef,type="l",xlab=expression(lambda),ylab="Coefficients",...)
	} else	if (xvar=="L1norm") {
		matplot(L1norm,coef,type="l",xlab=expression(sum(abs(beta[j]))),ylab="Coefficients",...)
		axis(4,at=coef[nrow(coef),],labels=paste(1:ncol(coef)))					
		} else {
			matplot(norm.frac,coef,type="l",xlab="|beta|/max|beta|",ylab="Coefficients",...)
			axis(4,at=coef[nrow(coef),],labels=paste(1:ncol(coef)))	
			}
	invisible()	
}
########################################################################################################################################################################
coef.nnlasso<-function(object,...)  object$coef
########################################################################################################################################################################
kfold<-function(data1,k)
{
	n=nrow(data1)
	p=ncol(data1)
	data1=cbind(rep(0,n),data1)
	#length of the folds
	l=floor(n/k)
	r=n%%k
	S=1:n
	for (i in 1:k)
	{
		if (i<k-r+1) 
		{
			s=sample(S,l,replace=F)
			data1[s,1]=i
		} else {
			s=sample(S,(l+1),replace=F)
			data1[s,1]=i
			}
		S=setdiff(S,s)		
	}
	return(data1)
}
########################################################################################################################################################################
fold<-function(data1,k,i)
{
	data1=kfold(data1,k)
	return(data1[which(data1[,1]==i),])
}
########################################################################################################################################################################
bars<-function(x, up, low, width = 0.03, ...) 
{
    xlim <- range(x)
    bw <- diff(xlim) * width
    segments(x, up, x, low, ...)
    segments(x - bw, up, x + bw, up, ...)
    segments(x - bw, low, x + bw, low, ...)
    range(up, low)
}
########################################################################################################################################################################
msefun.normal<-function(lambda1,f1,xi,yi)
{
	beta.new=predict.nnlasso(f1,mode="lambda",at=lambda1)	
	et=cbind(rep(1,nrow(xi)),xi)%*%beta.new
	return(sum((yi-et)^2)/length(yi))	
}
########################################################################################################################################################################
msefun.binomial<-function(lambda1,f1,xi,yi)
{
	beta.new=predict.nnlasso(f1,mode="lambda",at=lambda1)	
	et=cbind(rep(1,nrow(xi)),xi)%*%beta.new
	expet=exp(et)
	expet[which(expet==Inf)]=.Machine$double.xmax
	mu1=expet/(1+expet)
	t1=yi/mu1
	t2=(1-yi)/(1-mu1+.Machine$double.eps)
	dev=2*sum(yi*log(t1+.Machine$double.eps)+(1-yi)*log(t2+.Machine$double.eps))
	return(dev)	
}
########################################################################################################################################################################
msefun.poisson<-function(lambda1,f1,xi,yi)
{
	beta.new=predict.nnlasso(f1,mode="lambda",at=lambda1)	
	et=cbind(rep(1,nrow(xi)),xi)%*%beta.new
	mu1=exp(et)
	mu1[which(mu1==Inf)]=.Machine$double.xmax
	dev=2*sum(yi*log((yi+.Machine$double.eps)/(mu1+.Machine$double.eps))-(yi-mu1))
	return(dev)	
}
########################################################################################################################################################################
cv.nnlasso<-function(x,y,family=c("binomial","normal","poisson"),k=5,nlambda=50,tau=1,plot=TRUE,errorbars=TRUE)
{
  if(missing(family)) family="normal"
  switch(family, normal = cv.nnlasso.normal(x,y,k,nlambda=50,tau,plot,errorbars),binomial = cv.nnlasso.binomial(x,y,k,nlambda=50,tau,plot,errorbars), poisson= cv.nnlasso.poisson(x,y,k,nlambda=50,tau,plot,errorbars))
}
########################################################################################################################################################################
cv.nnlasso.normal<-function(x,y,k=5,nlambda=50,tau=1,plot=TRUE,errorbars=TRUE)
{
	meanx = colMeans(x)
        sx = scale(x, meanx, FALSE)
        meany = mean(y)
        sy = y - meany
    	one=rep(1,nrow(x))
	normx = sqrt(drop(one %*% (sx^2)))           
        sx = scale(sx, FALSE, normx)
    	mse=matrix(0,nlambda,k)
	max.lambda=max(abs(t(sx)%*%sy))
	rm(sx,sy)	
	gap=(log(1e-4)-log(max.lambda))/(nlambda-1)
	lambdas=exp(log(max.lambda)+gap*c(0:(nlambda-1)))
	x=kfold(x,k)
	folds.ident=x[,1]
	x=x[,2:ncol(x)]	
	for (i in 1:k)
	{
		ithfold=which(folds.ident==i)
		xi=x[ithfold,]
		yi=y[ithfold]
		x_i=x[-ithfold,]
		y_i=y[-ithfold]	
		f1=nnlasso.normal(x_i,y_i,tau=tau)
		mse[,i]=sapply(lambdas,msefun.normal,f1,xi,yi)				
	}
	se=sqrt(apply(mse,1,var)/k)
	pmse=rowMeans(mse)
	lambda=lambdas[which.min(pmse)]	
	out=list(lambda=lambda,pmse=pmse,lambdas=lambdas,se=se)
	if (plot) 
	{	
		plot(lambdas,pmse,xlab=expression(lambda),ylab="PMSE",ylim = range(pmse, pmse + se, pmse - se))
		if (errorbars) bars(lambdas, pmse+ se, pmse - se, width = 1/length(lambdas))
		abline(v=lambda)
	}
	return(out)
}	
########################################################################################################################################################################
cv.nnlasso.binomial<-function(x,y,k=5,nlambda=50,tau=1,plot=TRUE,errorbars=TRUE)
{
	one = rep(1, nrow(x))  	   	
    	meanx = colMeans(x)
        sx = scale(x, meanx, FALSE)        	
    	normx = sqrt(drop(one %*% (sx^2)))           
        sx = scale(sx, FALSE, normx)
    	max.lambda=max(abs(t(sx)%*%(y-0.5)))
	if (nrow(x)>ncol(x)) min.lambda=1e-4*max.lambda else min.lambda=1e-2*max.lambda	
	gap=(log(min.lambda)-log(max.lambda))/(nlambda-1)
	lambdas=exp(log(max.lambda)+gap*c(0:(nlambda-1)))
	mse=matrix(0,nlambda,k)
	x=kfold(x,k)
	folds.ident=x[,1]
	x=x[,2:ncol(x)]		
	for (i in 1:k)	
	{
		ithfold=which(folds.ident==i)
		xi=x[ithfold,]
		yi=y[ithfold]
		x_i=x[-ithfold,]
		y_i=y[-ithfold]			
		f1=nnlasso.binomial(x_i,y_i,tau=tau)
		mse[,i]=sapply(lambdas,msefun.binomial,f1,xi,yi)
	}
	se=sqrt(apply(mse,1,var)/k)
	pmse=rowMeans(mse)
	lambda=lambdas[which.min(pmse)]	
	out=list(lambda=lambda,pmse=pmse,lambdas=lambdas,se=se)
	if (plot) 
	{	
		plot(lambdas,pmse,xlab=expression(lambda),ylab="Deviance",ylim = range(pmse, pmse + se, pmse - se))
		if (errorbars) bars(lambdas, pmse+ se, pmse - se, width = 1/length(lambdas))
		abline(v=lambda)
	}
	return(out)
}
########################################################################################################################################################################
cv.nnlasso.poisson<-function(x,y,k=5,nlambda=50,tau=1,plot=TRUE,errorbars=TRUE)
{
	one = rep(1, nrow(x))  	   	
    	meanx = colMeans(x)
        sx = scale(x, meanx, FALSE)        	
    	normx = sqrt(drop(one %*% (sx^2)))           
        sx = scale(sx, FALSE, normx)
    	max.lambda=max(abs(t(sx)%*%(y-1)))
	if (nrow(x)>ncol(x)) min.lambda=1e-4*max.lambda else min.lambda=1e-2*max.lambda	
	gap=(log(min.lambda)-log(max.lambda))/(nlambda-1)
	lambdas=exp(log(max.lambda)+gap*c(0:(nlambda-1)))
	mse=matrix(0,nlambda,k)
	x=kfold(x,k)
	folds.ident=x[,1]
	x=x[,2:ncol(x)]		
	for (i in 1:k)	
	{
		ithfold=which(folds.ident==i)
		xi=x[ithfold,]
		yi=y[ithfold]
		x_i=x[-ithfold,]
		y_i=y[-ithfold]			
		f1=nnlasso.poisson(x_i,y_i,tau=tau)
		mse[,i]=sapply(lambdas,msefun.poisson,f1,xi,yi)
	}	
	se=sqrt(apply(mse,1,var)/k)
	pmse=rowMeans(mse)
	lambda=lambdas[which.min(pmse)]	
	out=list(lambda=lambda,pmse=pmse,lambdas=lambdas,se=se)
	if (plot) 
	{	
		plot(lambdas,pmse,xlab=expression(lambda),ylab="Deviance",ylim = range(pmse, pmse + se, pmse - se))
		if (errorbars) bars(lambdas, pmse+ se, pmse - se, width = 1/length(lambdas))
		abline(v=lambda)
	}
	return(out)
}	
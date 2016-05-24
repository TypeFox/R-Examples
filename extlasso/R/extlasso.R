extlasso<-function(x,y,family=c("normal","binomial","poisson"),intercept=TRUE,normalize=TRUE,tau=1,alpha=1e-12,eps=1e-6,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4)
{
	if(missing(family)) family="normal"	
  	switch(family, normal = extlasso.normal(x,y,intercept,normalize,tau,alpha,eps,tol,maxiter,nstep,min.lambda),binomial = extlasso.binomial(x,y,intercept,normalize,tau,alpha,eps,tol,maxiter,nstep,min.lambda), poisson = extlasso.poisson(x,y,intercept,normalize,tau,alpha,eps,tol,maxiter,nstep,min.lambda))
}
##############################################################################################################################################
extlasso.normal<-function(x,y,intercept=TRUE,normalize=TRUE,tau=1,alpha=1e-12,eps=1e-6,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4)
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
		dxpx=rep(1,p)
    	} else {
        	normx = rep(1, p)
		dxpx=colSums(x^2)
               }
	tx=t(x)
	xpy=tx%*%y
	if (n>p) xpx=tx%*%x else xpx=NULL
	max.lambda=max(abs(xpy))	
	stepsize=exp((log(min.lambda)-log(max.lambda))/nstep)
	lambdas=max.lambda*stepsize^((1:nstep)-1)
	coef=matrix(0,nstep,p)	
	of.value=rep(0,nstep)
	xbeta.old=rep(0,n)
	xpxbetaold=rep(0,p)
	of.value[1]=sum(y^2)/2
	lambda.iter=rep(0,nstep)
	g1=-xpy
	for(iter in 2:nstep)
	{
		active=which(coef[iter-1,] !=0 | abs(g1)>=tau*(2*lambdas[iter]-lambdas[iter-1]))	
		kkt=FALSE
		while(kkt==FALSE)
		{
			p1=length(active)
			res=extlasso.norm.lambda(n,p,p1,x[,active],y,xpx[active,active],dxpx[active],xpy[active],beta.old=coef[iter-1,active],tau,alpha,lambda1=lambdas[iter],tol,maxiter,eps,xbeta.old)
			if (res$conv=="yes")
			{	
				coef[iter,active]=res$beta.new
				if (n>p)
				{
					if (p1>1) xpxbetaold=xpx[,active]%*%res$beta.new else xpxbetaold=xpx[,active]*res$beta.new
				} else xpxbetaold=tx%*%res$xbeta.new
				g1=-xpy+xpxbetaold+lambdas[iter]*2*(1-tau)*coef[iter,]
				inactive=setdiff(1:p,active)
				indices=which(abs(g1[inactive])/(tau*lambdas[iter])>=0.99)
				if ((length(indices))==0)
				{
					xbeta.old=res$xbeta.new
					lambda.iter[iter]=res$iter
					of.value[iter]=res$ofv.new
					kkt=TRUE
				} else active=union(active,inactive[indices])										
			} else stop("The algorithm did not converge")
		}
	}
	coef=scale(coef,center=FALSE,scale=normx)
	if (intercept) beta0=rep(meany,nstep)-coef%*%meanx else beta0=rep(0,nstep)
	L1norm=rowSums(abs(coef))
	norm.frac=L1norm/max(L1norm)
	obj=list(beta0=beta0,coef=coef,lambdas=lambdas,L1norm=L1norm,norm.frac=norm.frac,lambda.iter=lambda.iter,of.value=of.value,normx=normx)
	class(obj)='extlasso'
	return(obj)
}
##############################################################################################################################################################
extlasso.binomial<-function(x,y,intercept=TRUE,normalize=TRUE,tau=1,alpha=1e-12,eps=1e-6,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4)
{
	np=dim(x)
	n=np[1]
	p=np[2]	 	
	if (intercept) 
	{
        	meanx = colMeans(x)
        	x = scale(x, meanx, FALSE)		
	} else meanx = rep(0, p)
    	if (normalize) 
  	{
        	normx = sqrt(colSums(x^2))           
        	x = scale(x, FALSE, normx)
	} else normx = rep(1, p)
	mtx=-t(x)
	tx2=mtx^2
	sumy=sum(y)
	beta0=rep(0,nstep)	
	coef=matrix(0,nstep,p)	
	of.value=rep(0,nstep)	
	xbeta.old=matrix(0,nstep,n)
	g1=mtx%*%(y-0.5)
	max.lambda=max(abs(g1))
	stepsize=exp((log(min.lambda)-log(max.lambda))/nstep)
	lambdas=max.lambda*stepsize^((1:nstep)-1)	
	lambda.iter=rep(0,nstep)				
	Mu=matrix(0.5,nstep,n)
	of.value[1]=n*log(2)
	for(iter in 2:nstep)
	{
		active=which(coef[iter-1,] !=0 | abs(g1)>=tau*(2*lambdas[iter]-lambdas[iter-1]))
		kkt=FALSE
		while(kkt==FALSE)
		{
			dxkx0=ifelse(intercept,sum(Mu[iter-1,]*(1-Mu[iter-1,])),Inf)
			dxkx1=tx2[active,]%*%(Mu[iter-1,]*(1-Mu[iter-1,]))
			p1=length(active)
			res=extlasso.binom.lambda(n,p,p1,sumy,beta0.old=beta0[iter-1],beta1.old=coef[iter-1,active],x[,active],y,dxkx0,dxkx1,tau,lambda1=lambdas[iter],alpha,tol,maxiter,eps,xbeta.old[iter-1,],mu1=Mu[iter-1,])
			if (res$conv=="yes") 
			{
				coef[iter,active]=res$beta1.new
				g1=mtx%*%(y-res$mu1)+2*lambdas[iter]*(1-tau)*coef[iter,]
				inactive=setdiff(1:p,active)
				indices=which(abs(g1[inactive])/(tau*lambdas[iter])>=0.99)
				if ((length(indices))==0)
				{	
					lambda.iter[iter]=res$iter
					of.value[iter]=res$ofv.new
					xbeta.old[iter,]=res$xbeta.new
					Mu[iter,]=res$mu1
					beta0[iter]=res$beta0.new
					kkt=TRUE
				} else active=union(active,inactive[indices])					
			} else stop("the algorithm did not converge")	
		}
	}	
	coef=scale(coef,center=FALSE,scale=normx)
	L1norm=rowSums(abs(coef))
	norm.frac=L1norm/max(L1norm)
	if (intercept) beta0=beta0-coef%*%meanx		
	obj=list(beta0=beta0,coef=coef,lambdas=lambdas,L1norm=L1norm,norm.frac=norm.frac,lambda.iter=lambda.iter,of.value=of.value,normx=normx)
	class(obj)='extlasso'
	return(obj)
}
##################################################################################################################################################################
extlasso.poisson<-function(x,y,intercept=TRUE,normalize=TRUE,tau=1,alpha=1e-12,eps=1e-6,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4)
{
	np=dim(x)
	n=np[1]
	p=np[2]	   	
    	if (intercept) 
	{
        	meanx = colMeans(x)
        	x = scale(x, meanx, FALSE)		   	
    	} else 	meanx = rep(0, p)
        if (normalize) 
  	{
        	normx = sqrt(colSums(x^2))           
        	x = scale(x, FALSE, normx)
	} else normx = rep(1, p)		
        mtx=-t(x)
	tx2=mtx^2
	sumy=sum(y)
	beta0=rep(0,nstep)   
	coef=matrix(0,nstep,p)	
	of.value=rep(0,nstep)	
	xbeta.old=matrix(0,nstep,n)
	g1=mtx%*%(y-1)	
	max.lambda=max(abs(g1))
	stepsize=exp((log(min.lambda)-log(max.lambda))/nstep)
	lambdas=max.lambda*stepsize^((1:nstep)-1)		
	lambda.iter=rep(0,nstep)				
	Mu=matrix(1,nstep,n)
	of.value[1]=n
	for(iter in 2:nstep)
	{		
		active=which(coef[iter-1,] !=0 | abs(g1)>=tau*(2*lambdas[iter]-lambdas[iter-1]))
		kkt=FALSE
		while(kkt==FALSE)
		{
			dxkx0=ifelse(intercept,sum(Mu[iter-1,]),Inf)
			dxkx1=tx2[active,]%*%Mu[iter-1,]
			p1=length(active)
			res=extlasso.pois.lambda(n,p,p1,sumy,beta0.old=beta0[iter-1],beta1.old=coef[iter-1,active],x[,active],y,dxkx0,dxkx1,tau,lambda1=lambdas[iter],alpha,tol,maxiter,eps,xbeta.old[iter-1,],mu1=Mu[iter-1,])
			if (res$conv=="yes") 
			{
				coef[iter,active]=res$beta1.new
				g1=mtx%*%(y-res$mu1)+2*lambdas[iter]*(1-tau)*coef[iter,] 
				inactive=setdiff(1:p,active)
				indices=which(abs(g1[inactive])/(tau*lambdas[iter])>=0.99)
				if ((length(indices))==0)
				{
					lambda.iter[iter]=res$iter
					of.value[iter]=res$ofv.new
					xbeta.old[iter,]=res$xbeta.new
					Mu[iter,]=res$mu1
					beta0[iter]=res$beta0.new
					kkt=TRUE
				} else active=union(active,inactive[indices])
			} else stop("the algorithm did not converge")	
		}
	}
	coef=scale(coef,center=FALSE,scale=normx)
	L1norm=rowSums(abs(coef))
	norm.frac=L1norm/max(L1norm)
	if (intercept) beta0=beta0-coef%*%meanx		
	obj=list(beta0=beta0,coef=coef,lambdas=lambdas,L1norm=L1norm,norm.frac=norm.frac,lambda.iter=lambda.iter,of.value=of.value,normx=normx)
	class(obj)='extlasso'
	return(obj)
}
#################################################################################################################################################################
extlasso.norm.lambda<-function(n,p,p1,x,y,xpx,dxpx,xpy,beta.old,tau,alpha,lambda1,tol,maxiter,eps,xbeta.old)
{
	epp=0.001
	if (n<=p) tx=t(x)	
	ofv.old=sum((y-xbeta.old)^2)/2+lambda1*(tau*(sum(sqrt(beta.old^2+alpha)))+(1-tau)*sum(beta.old^2))
	for (iter in 1:maxiter)
	{
		if (n>p)
		{
			if (p1>1) xpxbetaold=xpx%*%beta.old else xpxbetaold=xpx*beta.old
		} else xpxbetaold=tx%*%xbeta.old
		apbeta=sqrt(beta.old^2+alpha)		
		g=-xpy+xpxbetaold+lambda1*(tau*beta.old/apbeta+2*(1-tau)*beta.old) 
		b=g/(dxpx+lambda1*(tau*alpha/(apbeta^3)+2*(1-tau)*rep(1,p1)))
		beta.new=beta.old-b
		if (p1>1) xb=x%*%b else xb=x*b
		xbeta.new=xbeta.old-xb					       
		ofv.new=sum((y-xbeta.new)^2)/2+lambda1*(tau*(sum(sqrt(beta.new^2+alpha)))+(1-tau)*sum(beta.new^2))
		delta=1
		t1=epp*sum(g*b)		
		while (ofv.new-delta*t1>ofv.old & delta>1e-5)			
		{	
			delta=delta/2
			beta.new=beta.old-delta*b			
			xbeta.new=xbeta.old-delta*xb
			ofv.new=sum((y-xbeta.new)^2)/2+lambda1*(tau*(sum(sqrt(beta.new^2+alpha)))+(1-tau)*sum(beta.new^2))
		}
		if (ofv.new-delta*t1>ofv.old & delta<=1e-5) 
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
	if (iter<maxiter)
	{
		conv="yes"		
		beta.new[which(abs(beta.new)<eps)]=0		
	} else conv="no"
	res=list(beta.new=beta.new,conv=conv,iter=iter,ofv.new=ofv.new,xbeta.new=xbeta.new) 
	return(res)		
}
####################################################################################################################################################
extlasso.binom.lambda<-function(n,p,p1,sumy,beta0.old,beta1.old,x,y,dxkx0,dxkx1,tau,lambda1,alpha,tol,maxiter,eps,xbeta.old,mu1)
{
	epp=0.001
	mtx=-t(x)
	expeta.old=exp(xbeta.old)
	ofv.old=-(sum(y*xbeta.old)-sum(log(1+expeta.old)))+lambda1*(tau*sum(sqrt(beta1.old^2+alpha))+(1-tau)*sum(beta1.old^2))
	for (iter in 1:maxiter)
	{
		g0=sum(mu1)-sumy
		b0=g0/dxkx0
		beta0.new=beta0.old-b0
		apbeta=sqrt(beta1.old^2+alpha)
		g1=mtx%*%(y-mu1)+lambda1*(tau*beta1.old/apbeta+2*(1-tau)*beta1.old)
		b1=g1/(dxkx1+lambda1*(tau*alpha/(apbeta^3)+2*(1-tau)))
		beta1.new=beta1.old-b1
		t1=epp*(sum(g1*b1)+g0*b0)
		if (p1>1) xb=x%*%b1+b0 else xb=x*b1+b0		
		xbeta.new=xbeta.old-xb
		expeta.new=exp(xbeta.new)
		ofv.new=-(sum(y*xbeta.new)-sum(log(1+expeta.new)))+lambda1*(tau*sum(sqrt(beta1.new^2+alpha))+(1-tau)*sum(beta1.new^2))
		delta=1
		while (ofv.new-delta*t1>ofv.old & delta>1e-5)			
		{	
			delta=delta/2
			beta0.new=beta0.old-delta*b0 
			beta1.new=beta1.old-delta*b1			 
			xbeta.new=xbeta.old-delta*xb
			expeta.new=exp(xbeta.new)
			ofv.new=-(sum(y*xbeta.new)-sum(log(1+expeta.new)))+lambda1*(tau*sum(sqrt(beta1.new^2+alpha))+(1-tau)*sum(beta1.new^2))
		}
		if (ofv.new-delta*t1>ofv.old & delta<1e-5) 
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
	if (iter<maxiter) 
	{
		conv="yes"		
		beta1.new[which(abs(beta1.new)<eps)]=0
	} else conv="no"	
	res=list(beta0.new=beta0.new,beta1.new=beta1.new,conv=conv,iter=iter,ofv.new=ofv.new,xbeta.new=xbeta.new,mu1=mu1) 
	return(res)
}
##############################################################################################################################################################
extlasso.pois.lambda<-function(n,p,p1,sumy,beta0.old,beta1.old,x,y,dxkx0,dxkx1,tau,lambda1,alpha,tol,maxiter,eps,xbeta.old,mu1)
{
	epp=0.001
	mtx=-t(x)
	expeta.old=exp(xbeta.old)
	ofv.old=-(sum(y*xbeta.old)-sum(expeta.old))+lambda1*(tau*sum(sqrt(beta1.old^2+alpha))+(1-tau)*sum(beta1.old^2))	
	for (iter in 1:maxiter)
	{
		g0=sum(mu1)-sumy
		b0=g0/dxkx0
		beta0.new=beta0.old-b0
		apbeta=sqrt(beta1.old^2+alpha)
		g1=mtx%*%(y-mu1)+lambda1*(tau*beta1.old/apbeta+2*(1-tau)*beta1.old)
		b1=g1/(dxkx1+lambda1*(tau*alpha/(apbeta^3)+2*(1-tau)))
		beta1.new=beta1.old-b1
		t1=epp*(sum(g1*b1)+g0*b0)
		if (p1>1) xb=x%*%b1+b0 else xb=x*b1+b0
		xbeta.new=xbeta.old-xb
		expeta.new=exp(xbeta.new)
		ofv.new=-(sum(y*xbeta.new)-sum(expeta.new))+lambda1*(tau*sum(sqrt(beta1.new^2+alpha))+(1-tau)*sum(beta1.new^2))	
		delta=1		
		while (ofv.new-delta*t1>ofv.old & delta>1e-5)			
		{	
			delta=delta/2
			beta0.new=beta0.old-delta*b0
			beta1.new=beta1.old-delta*b1		
			xbeta.new=xbeta.old-delta*xb
			expeta.new=exp(xbeta.new)
			ofv.new=-(sum(y*xbeta.new)-sum(expeta.new))+lambda1*(tau*(sum(sqrt(beta1.new^2+alpha)))+(1-tau)*sum(beta1.new^2))	
		}
		if (ofv.new-delta*t1>ofv.old & delta<1e-5) 
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
	if (iter<maxiter) 
	{
		conv="yes"		
		beta1.new[which(abs(beta1.new)<eps)]=0
	} else conv="no"	
	res=list(beta0.new=beta0.new,beta1.new=beta1.new,conv=conv,iter=iter,ofv.new=ofv.new,xbeta.new=xbeta.new,mu1=mu1) 
	return(res)
}
###############################################################################################################################################################################
predict.extlasso<-function(object,mode=c("fraction","norm","lambda"),at=0,...)
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
coef.extlasso<-function(object,...)  object$coef
########################################################################################################################################################################
plot.extlasso<-function(x,xvar=c("lambda","L1norm","fraction of norm"),...)
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
	beta.new=predict.extlasso(f1,mode="lambda",at=lambda1)	
	et=cbind(rep(1,nrow(xi)),xi)%*%beta.new
	return(sum((yi-et)^2)/length(yi))	
}
########################################################################################################################################################################
msefun.binomial<-function(lambda1,f1,xi,yi)
{
	beta.new=predict.extlasso(f1,mode="lambda",at=lambda1)	
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
	beta.new=predict.extlasso(f1,mode="lambda",at=lambda1)	
	et=cbind(rep(1,nrow(xi)),xi)%*%beta.new
	mu1=exp(et)
	mu1[which(mu1==Inf)]=.Machine$double.xmax
	dev=2*sum(yi*log((yi+.Machine$double.eps)/(mu1+.Machine$double.eps))-(yi-mu1))
	return(dev)	
}
########################################################################################################################################################################
cv.extlasso<-function(x,y,family=c("binomial","normal","poisson"),k=5,nlambda=50,tau=1,plot=TRUE,errorbars=TRUE)
{
  if(missing(family)) family="normal"
  switch(family, normal = cv.normal(x,y,k,nlambda=50,tau,plot,errorbars),binomial = cv.binomial(x,y,k,nlambda=50,tau,plot,errorbars), poisson= cv.poisson(x,y,k,nlambda=50,tau,plot,errorbars))
}
########################################################################################################################################################################
cv.normal<-function(x,y,k=5,nlambda=50,tau=1,plot=TRUE,errorbars=TRUE)
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
	gap=exp((log(1e-4)-log(max.lambda))/nlambda)
	lambdas=max.lambda*gap^((1:nlambda)-1)
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
		f1=extlasso.normal(x_i,y_i,tau=tau)
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
cv.binomial<-function(x,y,k=5,nlambda=50,tau=1,plot=TRUE,errorbars=TRUE)
{
	one = rep(1, nrow(x))  	   	
  	meanx = colMeans(x)
	sx = scale(x, meanx, FALSE)        	
	normx = sqrt(drop(one %*% (sx^2)))           
  	sx = scale(sx, FALSE, normx)
  	max.lambda=max(abs(t(sx)%*%(y-0.5)))
	if (nrow(x)>ncol(x)) min.lambda=1e-4*max.lambda else min.lambda=1e-2*max.lambda
	gap=exp((log(min.lambda)-log(max.lambda))/nlambda)
	lambdas=max.lambda*gap^((1:nlambda)-1)
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
	  f1=extlasso.binomial(x_i,y_i,tau=tau)
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
cv.poisson<-function(x,y,k=5,nlambda=50,tau=1,plot=TRUE,errorbars=TRUE)
{
	one = rep(1, nrow(x))  	   	
    	meanx = colMeans(x)
        sx = scale(x, meanx, FALSE)        	
    	normx = sqrt(drop(one %*% (sx^2)))           
        sx = scale(sx, FALSE, normx)
    	max.lambda=max(abs(t(sx)%*%(y-1)))
	if (nrow(x)>ncol(x)) min.lambda=1e-4*max.lambda else min.lambda=1e-2*max.lambda
	gap=exp((log(min.lambda)-log(max.lambda))/nlambda)
	lambdas=max.lambda*gap^((1:nlambda)-1)
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
		f1=extlasso.poisson(x_i,y_i,tau=tau)
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
########################################################################################################################################################################
#extlasso <- function(x, ...) UseMethod("extlasso")
########################################################################################################################################################################
fl.lambda=function(n,p,x,y,xpx,dxpx,xpy,beta.old,ofv.old,alpha,lambda1,lambda2,tol,maxiter,eps,xbeta.old)
{
	epp=0.00001
	if (n<=p) tx=t(x)
	temp=rep(0,(p+1))
	betadiff.old=beta.old[2:p]-beta.old[1:(p-1)]
	for (iter in 1:maxiter)
	{
		if (n>p) xpxbetaold=xpx%*%beta.old else xpxbetaold=tx%*%xbeta.old
		apbeta=sqrt(beta.old^2+alpha)
		apbetadiff=sqrt(betadiff.old^2+alpha)
		temp[2:p]=betadiff.old/apbetadiff
		g=-xpy+xpxbetaold+lambda1*beta.old/apbeta+lambda2*(temp[1:p]-temp[2:(p+1)])
		temp[2:p]=alpha/(apbetadiff^3)
		b=g/(dxpx+lambda1*alpha/(apbeta^3)+lambda2*(temp[1:p]+temp[2:(p+1)]))
		beta.new=beta.old-b
		betadiff.new=beta.new[2:p]-beta.new[1:(p-1)]
		xb=x%*%b
		xbeta.new=xbeta.old-xb
		ofv.new=sum((y-xbeta.new)^2)/2+lambda1*(sum(sqrt(beta.new^2+alpha)))+lambda2*(sum(sqrt(betadiff.new^2+alpha)))
		delta=1
		t1=epp*(sum(g*b))		
		while (ofv.new-delta*t1>ofv.old & delta>1e-5)			
		{	
			delta=delta/2
			beta.new=beta.old-delta*b
			betadiff.new=beta.new[2:p]-beta.new[1:(p-1)]			
			xbeta.new=xbeta.old-delta*xb
			ofv.new=sum((y-xbeta.new)^2)/2+lambda1*(sum(sqrt(beta.new^2+alpha)))+lambda2*(sum(sqrt(betadiff.new^2+alpha)))
		}
		if (ofv.new-delta*t1>ofv.old & delta<=1e-5) 
		{
			beta.new=beta.old
			ofv.new=ofv.old	
			break
		}
		if(abs(ofv.old-ofv.new)<=tol) break
		beta.old=beta.new
		betadiff.old=betadiff.new		
		xbeta.old=xbeta.new
		ofv.old=ofv.new		
	}
	if (iter<maxiter)
	{
		conv="yes"		
		beta.new[which(abs(beta.new)<eps)]=0		
	} else conv="no"
	res=list(beta.new=beta.new,conv=conv,iter=iter,ofv.new=ofv.new) 
	return(res)		
}
########################################################################################################################################################################
fusedlasso=function(x,y,lambda1,lambda2,intercept=TRUE,normalize=TRUE,alpha=1e-6,eps=1e-6,tol=1e-8,maxiter=1e5)
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
		dxpx=rep(1,p)
    	} else {
        	normx = rep(1, p)
		dxpx=colSums(x^2)
              }
	tx=t(x)
	xpy=tx%*%y
	if (n>p) xpx=tx%*%x else xpx=NULL
	beta.old=rep(0,p)
	xbeta.old=rep(0,n)
	ofv.old=sum(y^2)/2
	res=fl.lambda(n,p,x,y,xpx,dxpx,xpy,beta.old,ofv.old,alpha,lambda1,lambda2,tol,maxiter,eps,xbeta.old)
	if (res$conv=="yes")
	{	coef=res$beta.new			
		lambda.iter=res$iter
		of.value=res$ofv.new			
	} else stop("The algorithm did not converge")
	dim(coef)=c(1,p)	
	coef=scale(coef,center=FALSE,scale=normx)
	L1norm=rowSums(abs(coef))
	if(intercept) beta0=meany-coef%*%meanx else beta0=0
	obj=list(beta0=beta0,coef=coef,lambda1=lambda1,lambda2=lambda2,L1norm=L1norm,lambda.iter=lambda.iter,of.value=of.value)
	class(obj)='extlasso'
	return(obj)
}
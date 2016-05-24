rm(list=ls())
library(tmvtnorm)
library(mvtnorm)
library(MCMCpack)
library(MASS)


sumwt <- function(x=NULL, z=NULL ) { sum(x*z) }

Zig <- function(x=NULL, gpar=NULL, v=NULL) {
	n = nrow(x)
	g = length(gpar$pi);
	if (g > 1) {
		zlog = glogskewtden(x=x, gpar=gpar) #I think that v is for annealing

		w = t(apply( zlog, 1, function(z,wt,v=1) { 
			x= exp( v*(z + log(wt)) ) 
			x=x/sum(x);
			return(x) }, wt=gpar$pi,v=v ))
	} else w = matrix(1,nrow=n, ncol=g)
	return(w)
	}

zrand <- function(x=NULL, G=NULL)	{
	priors=rep(1,G)
	z=rdirichlet(nrow(x), priors)
	return(z)
}

gparkmeans <- function(x=NULL, G=NULL){
	
	class.vec <- kmeans(x,G)$cluster
	z <- matrix(0, nrow=nrow(x), ncol=G)
	for(g in 1:G) z[,g] = as.numeric(class.vec==g) 
	return(z)
}

gparrand <- function(x=NULL, G=NULL){
	
	priors=rep(1,G)
	z=rdirichlet(nrow(x), priors)
	return(z)
}

rpar <- function(x=NULL, p=NULL, q=NULL, dlam =FALSE, z=NULL) {
	val = list()
	temp = cov.wt(x, wt=z)
	val$xi  = temp$center
	if ( dlam ) val$lam = diag(rnorm(p, 0, sd=.01))
	else val$lam = matrix(rnorm(p, 0, sd=.01),nrow=p, ncol=p) 
	val$sig = temp$cov

	val$invSig = solve(val$sig)
	p = nrow(val$lam)
	val$omega  = val$sig + t(val$lam) %*% val$lam
	val$invOmg = solve( val$omega )
	val$delta  = diag(p) - val$lam %*% val$invOmg %*% t(val$lam)
	
	esigma = eigen(val$sig)
	if(q==1) val$Lambda = t(sweep(t(esigma$vectors[,1:q]), 2, sqrt(esigma$values[1:q]), FUN="*"))
	if(q!=1) val$Lambda = sweep(esigma$vectors[,1:q], 2, sqrt(esigma$values[1:q]), FUN="*")
	val$psi = diag(c(diag(val$sig -  val$Lambda %*% t(val$Lambda))))
	if(p>1) P=solve(val$psi)
	if(p==1) P=1/val$psi
	val$invSig = P-P%*%val$Lambda%*%(diag(q)+t(val$Lambda)%*%val$psi%*%val$Lambda)%*%t(val$Lambda)%*%P


	val$nu   = 10
	return(val)
}

rgpar <- function(x=NULL, g=NULL, p=NULL, dlam=FALSE, q=NULL, z=NULL) {
	val = list()
	for (i in 1:g) val[[i]] = rpar(x=x, p=p, q=q, dlam = dlam, z=z[,i]) 
	val$pi = apply(z,2,mean)
	return(val)
}

aitken<-function(l, it){
	flag=1
	if(l[it-1] == l[it-2]){
		flag<-0
	}
	else{
		ak<-(l[it]-l[it-1])/(l[it-1]-l[it-2])
		l.inf<-l[it-1]+((l[it]-l[it-1])/(1-ak))
		if( abs(l.inf-l[it])< 0.1 ) flag<-0
		else flag<-1
	}	
	return(flag)
}

bic <- function(l=NULL, par=NULL, x=NULL, q=NULL){

	it=length(l)
	l.max=l[it]
	N=nrow(x)
	G=length(par)
	p=ncol(x)
	m=G*(3*p+p*q+1)

	val = 2*l.max - m*log(N)
	
	return(val)
}

uskewFA <- function(x=NULL, G=NULL, q=NULL, init=1, max.it=100) {
	
	x=as.matrix(x)
	val = list()
	val$z = matrix(0,nrow=nrow(x),ncol=G)
	if(init==1) val$z=gparkmeans(x=x, G=G)
	if(init==2) val$z=gparrand(x=x, G=G) 
	
	val$gpar = rgpar(x=x, g=G, p=ncol(x), dlam=TRUE, q=q, z=val$z)
	n=nrow(x)
	it=1
	not.converged=1 
	while(not.converged && it<=100) {		
		val$z = Zig(x=x, gpar=val$gpar, v=1)
		for (g in 1:G) {
 		estepg = eStep(x=x, par=val$gpar[[g]], wt=val$z[,g])
 		val$gpar[[g]] = mStep(x= x, gpar=val$gpar[[g]], estep=estepg, z=val$z[,g], it=it)
 		}
 		val$gpar$pi = apply(val$z,2,mean)
 		#print(val$gpar$pi)
 		
		val$loglik[it] = loglik(x=x, gpar=val$gpar)
		if(it>3) not.converged = aitken(l=val$loglik, it=it)
		#print(it)
		#plot(val$loglik)
		it=it+1
	}
	val$z   = Zig(x=x, gpar=val$gpar, v=1)
	val$map = apply(val$z,1, function(z){ (1:length(z))[z==max(z)] }) 
	output=list()
	output$map=val$map
	output$bic=bic(l=val$loglik,par=val$gpar,x=x,q=q)
	output$zhat=val$z
	output$likelihood=val$loglik
	return(val)
	}
	
loglik <- function(x=NULL, gpar=NULL) {
	# output is a G x nrow(data) matrix
	zlog = glogskewtden(x= x, gpar= gpar)
	
	w = apply( exp(zlog),1,function(z,wt) { sum(z*wt) } , wt=gpar$pi)
	val = sum(log(w))
	if( is.nan(val) ) val =NA	
	return(val)
	}

#Returns log of density to be used to calculate log likelihood
glogskewtden <- function(x=NULL, gpar=NULL) {
	n = nrow(x); g = length(gpar$pi);
	zlog = matrix(0, nrow=n, ncol=g)
	for (k in 1:g) zlog[,k] = skewtden(x=x, par= gpar[[k]], logd= TRUE)
	return(zlog)	
}

updateNu <- function(nuv=NULL) {
	nu0 = c(2,200)
	v0 = log(nu0/2) + 1 - digamma(nu0/2) - nuv
	sv0 = sign(v0)
	
	if ( all(sv0<0) ) { 
		val = 2
	} else if ( all(sv0>0) ) {
		val= 200
	} else {
		temp = uniroot(f=function(nu=NULL, nuv=NULL) { 
		val = log(nu/2) + 1 - digamma(nu/2) - nuv
		return(val) }, interval=c(2,200), nuv=nuv)
		val = temp$root
	} 
	return(val)
}

eStep <- function(x=NULL, par=NULL,wt=NULL) {
	n = nrow(x); p =ncol(x); 
	if (is.null(wt)) wt= rep(1,n)
	
	eta = matrix(0, nrow=n, ncol=p)
	psi = matrix(0, nrow=p, ncol=p)
	tau = numeric(n)
	logtau = numeric(n)
	u = mahalanobis(x, center=par$xi, cov= par$invOmg, inverted=TRUE)
	qq = as.matrix(sweep(x, 2, par$xi, "-")) %*% par$invOmg %*% t(par$lam)
	lx= rep(-Inf, p)

	for (i in 1:nrow(x)) {
		cc = (par$nu + p + c(0,2))/(par$nu + u[i])
		pp = c(	pmvt2(lower=lx, upper=qq[i,]*sqrt(cc[1]), delta=rep(0,p), sigma=round(par$delta,3), df=par$nu+p),
				pmvt2(lower=lx, upper=qq[i,]*sqrt(cc[2]), delta=rep(0,p), sigma=round(par$delta,3), df=par$nu+p+2) ) 
		tau[i] = cc[1] * pp[2]/pp[1]

# if (tau[i] <0) {
# print( c(tau[i], cc,pp) )
	# print('tau <0')
	# print(lx)
	# print(qq[i,])
	# print(par$delta)
	# print(par$nu+p+2)
	# #print(par)
	# stop('here')	
# }		
		logtau[i] = tau[i] - log( (par$nu+u[i])/2) - cc[1] + digamma( (par$nu+p)/2 )
		
		
		temp=try(truncatedTmom(mu=qq[i,], sigma=round(par$delta/cc[2]), a=rep(0,p), nu=par$nu + p + 2 ),silent=TRUE)
		if(is.list(temp)) temp=temp else temp = try(truncatedTmom(mu=qq[i,], sigma=round(par$delta/cc[2],3), a=rep(0,p), nu=par$nu + p + 2 ),silent=TRUE)
		if(is.list(temp)) temp=temp else temp = try(truncatedTmom(mu=qq[i,], sigma=round(par$delta/cc[2],2), a=rep(0,p), nu=par$nu + p + 2 ),silent=TRUE)
		
		eta[i,]  = temp$tmean
		psi      = psi + wt[i]*tau[i]*( temp$tvar + outer(temp$tmean,temp$tmean) )
	}
	
	val = list(eta=eta, psi.wt=psi, tau=tau, logtau=logtau)
	return(val)
}


mStep <- function(x= NULL, gpar=NULL, estep=NULL, z=NULL, it=NULL) {
	n = nrow(x); p=ncol(x); q=ncol(gpar$Lambda); #

	wtt = z*estep$tau
	ngt = sum(wtt)
	par1=gpar

	par1$xi = as.numeric(apply(x - estep$eta %*% par1$lam,2, weighted.mean, w=wtt))
	r = sweep(x,2, par1$xi, FUN="-") 
	r = as.matrix(r,nrow=n,ncol=p,dimnames=NULL)
	# diagonal lam
	temp = matrix(0,nrow=p,ncol=p)
	for (j in 1:n) temp = temp+ wtt[j]*outer(as.vector(r[j,]), as.vector(estep$eta[j,]))
	
# # 	 if(!is.null(par1$Delta)) par1$lam = diag( as.numeric(solve(par1$invSig * estep$psi.wt) %*%  (par1$invSig * t(temp))  %*% matrix(1,nrow=p,ncol=1)), p, p)	
	 # else par1$Delta = solve(estep$psi.wt) %*% t(temp) 
	 par1$lam = par1$lam = diag( as.numeric(solve(par1$invSig * estep$psi.wt) %*%  (par1$invSig * t(temp))  %*% matrix(1,nrow=p,ncol=1)), p, p)

	rr = cov.wt(r, wt=wtt, center=rep(0,p), method="ML")$cov* ngt
	sig = rr + t(par1$lam) %*% (estep$psi.wt) %*% par1$lam  - ( (temp) %*% (par1$lam)  + t(par1$lam) %*% t(temp) )
	par1$sig = sig/sum(wtt)

	par1$nu = updateNu(nuv=weighted.mean( estep$tau - estep$logtau , w=z) )
	par1$Beta = t(par1$Lambda)%*%par1$invSig
	par1$omega  = par1$sig + t(par1$lam) %*% par1$lam;
	par1$invOmg = solve( par1$omega ); 
	par1$psi = diag(diag(par1$sig-par1$Lambda%*%par1$Beta%*%par1$sig))
	#par1$inv.psi = solve(par1$psi)
	old.Theta=par1$Theta
	par1$Theta = diag(q)-par1$Beta%*%par1$Lambda+par1$Beta%*%par1$sig%*%t(par1$Beta)
	#if(any(par1$Theta>10) && it>=4) par1$Theta=old.Theta
	#print(par1$Theta)
	par1$Lambda = par1$sig%*%t(par1$Beta)%*%ginv(par1$Theta)
	
	#P=par1$inv.psi
	par1$invSig = par1$invOmg #P-P%*%par1$Lambda%*%(diag(q)+t(par1$Lambda)%*%par1$psi%*%par1$Lambda)%*%t(par1$Lambda)%*%P
	
	#par1$delta  = diag(p) - (par1$lam) %*% par1$invOmg %*% t(par1$lam);
	return(par1)	
}

truncatedTmom <- function(mu=NULL, sigma=NULL, a=NULL, nu=NULL) {
	if (length(mu) == 1 ) val = truncatedTmom1(mu= as.numeric(mu), sigma= as.numeric(sigma), lower=0, upper=Inf, nu=nu)
	else if (length(mu) == 2) val = truncatedTmom2(mu= mu, sigma= sigma, a=a, nu=nu)
	else val = truncatedTmomk(mu= mu, sigma= sigma, a=a, nu=nu)
	return(val)
}

truncatedTmomk <- function(mu=NULL, sigma=NULL, a=NULL, nu=NULL) {
	mu = -1*mu
	p = length(mu)
	if (p < 3 ) stop("truncatedTmomk is for k greater 2")
		
	xi = numeric(p)
	H  = matrix(0,p,p)
	for (k in 1:p) {
		lv1 = 1/2*( log(nu/2) - log(2*pi*sigma[k,k]) )
		lv2 = (nu-1)/2*( log(nu)-log(nu+(mu[k]-a[k])^2/sigma[k,k]) )
		lv3 = lgamma((nu-1)/2) - lgamma(nu/2)
		
		a.star = (a[-k]-mu[-k]) -  ((a[k]-mu[k])/sigma[k,k]) * sigma[k,-k] 
		s.star = (nu + (mu[k]-a[k])^2/sigma[k,k])/(nu-1) *( sigma[-k,-k] - outer(sigma[k,-k],sigma[k,-k])/sigma[k,k] )
		
		v4 = pmvt2(lower=rep(-Inf,p-1), upper=a.star, delta=rep(0,p-1), sigma=s.star, df=nu-1)
		xi[k] = exp(lv1+lv2+lv3) *v4
		
		for (l in 1:p) {
		if (k != l ) {
			akl  = a[c(k,l)];  a.kl  = a[-c(k,l)];
			mukl = mu[c(k,l)]; mu.kl = mu[-c(k,l)];
			slk  = sigma[c(k,l),c(k,l)]; s.lk  = sigma[-c(k,l),-c(k,l)];
			s.lklk	= as.matrix(sigma[c(k,l), -c(k,l)])
						
			invslk = solve(slk)
			
			nu.star = nu + as.numeric( (akl - mukl) %*% invslk %*% (akl - mukl) )
			a.ss = as.numeric( (a.kl - mu.kl) - t(s.lklk) %*% invslk %*% (akl - mukl) )
			s.ss = nu.star/(nu-2)*( s.lk - t(s.lklk) %*% invslk %*% (s.lklk) ) 

			lv1 = log(nu/(nu-2))-log(2*pi)-log( sigma[k,k]*sigma[l,l] - sigma[k,l]*sigma[l,k] )/2
			lv2 = (nu/2-1)*log(nu/nu.star)
			v3  = pmvt2(lower=rep(-Inf,p-2), upper=a.ss, delta=rep(0,p-2), sigma=s.ss, df=nu-2)

			H[k,l] = -exp(lv1+lv2)*v3
		}} 		
	}		
	for (k in 1:p) H[k,k] = ( xi[k]*(a[k]-mu[k]) - sum(sigma[-k,k]*H[-k,k]))/sigma[k,k]	
		
	c1 = pmvt2(lower=rep(-Inf,p), upper=a-mu, delta=rep(0,p), sigma=sigma, df=nu)
	c2 = pmvt2(lower=rep(-Inf,p), upper=a-mu, delta=rep(0,p), sigma=sigma*(nu/(nu-2)), df=nu-2)
	
	mu.star = as.numeric( (sigma %*% xi)/c1 )
	tmu = mu - mu.star
		
	s2 = ( sigma %*% H %*% sigma )/c1
	s3 = ( c2/c1*nu/(nu-2) )* sigma
	tsig = s3 - s2 
	tsig = tsig - outer(mu-tmu, mu-tmu)

	
	val = list(tmean=-1*tmu, tvar=tsig)
	return(val)
}


truncatedTmom1 <- function(mu=NULL, sigma=NULL, lower=NULL, upper=NULL, nu=NULL) {
	s = sqrt(sigma)
	
	ab = (c(lower,upper) - mu)/s 

	v     = nu
	d1    = -dt( ab*sqrt((v-2)/v), df=nu-2)*sqrt(v/(v-2))
	p0    = pt( ab, df=nu)
	rdp   = diff(d1)/diff(p0) 
	mus   = s*rdp
	tmean = mu + mus

	p2  = pt( ab/sqrt( nu/(nu-2) ), df=nu-2)
	tvar  = ((nu-1)*(nu/(nu-2))*diff(p2)/diff(p0) - nu )*sigma
	#tvar = tvar  - (mu^2 - tmean*mu - tmean* mu ) - tmean^2
	tvar = tvar  - (mu- tmean)^2
	
	val = list(tmean=as.numeric(tmean), tvar=as.numeric(tvar)  )
	return(val)
}


truncatedTmom2 <- function(mu=NULL, sigma=NULL, a=NULL, nu=NULL) {
	mu = -1*mu
	
	if (length(mu) != 2 ) stop("mu does not have length 2")
	p = length(mu)
		
	xi = numeric(2)
	H  = matrix(0,2,2)
	
	ds = diag(sigma)
	cs = sigma[1,2]
		
	lv1 = 1/2*( log(nu/2) - log(2*pi*ds) )
	lv2 = (nu-1)/2*( log(nu)-log(nu+(mu-a)^2/ds) )
	lv3 = lgamma( (nu-1)/2) - lgamma(nu/2)

	a.star = rev(a-mu) -  ((a-mu)/ds)*cs
	s.star = ( rev(ds) - cs^2/ds )*(nu + (mu-a)^2/ds )/(nu-1) 
	v4 = pt(a.star/sqrt(s.star), df=nu-1)
	xi = exp(lv1+lv2+lv3) *v4

	nu.star = nu + as.numeric( (a - mu) %*% solve(sigma) %*% (a - mu) )
	
	lv1 = log(nu/(nu-2))-log(2*pi)-log( prod(ds) - cs^2 )/2
	lv2 = (nu/2-1)*log(nu/nu.star)
	H[1,2] = -exp(lv1+lv2)
	H[2,1] = H[1,2]

	diag(H) = (xi*(a-mu) - cs*H[1,2])/ds
##########
	c1 = pmvt2(lower=rep(-Inf,2), upper=a-mu, delta=rep(0,2), sigma=sigma, df=nu)
	c2 = pmvt2(lower=rep(-Inf,2), upper=a-mu, delta=rep(0,2), sigma=sigma*(nu/(nu-2)), df=nu-2)
	
	mu.star = as.numeric( (sigma %*% xi)/c1 )
	tmu = mu - mu.star
		
#	s1 = outer(mu-mu.star, mu-mu.star) - outer(mu.star,mu.star) 
#	s1 = outer(mu,mu)-outer(mu,mu.star)-outer(mu.star,mu)
	s2 = ( sigma %*% H %*% sigma )/c1
	s3 = ( c2/c1*nu/(nu-2) )* sigma

	tsig = s3 - s2 
	tsig = tsig - outer(mu-tmu, mu-tmu)
	
	val = list(tmean=-1*tmu, tvar=tsig)
	return(val)
}

skewtden <- function(x=NULL, par=NULL, logd=TRUE) {
	## x is the data
	n = nrow(x); p = ncol(x); 
	
	r  = sweep(x, 2, par$xi, "-")
	qq = as.matrix(r) %*% par$invOmg %*% t(par$lam)
	u = mahalanobis(r, center=rep(0,p), cov= par$invOmg, inverted=TRUE)

	u.star = sqrt( (p + par$nu)/(u + par$nu) )
	q.star = sweep(qq, 1, u.star, FUN="*")

	v1 = p*log(2)
	v2= try(dmvt(r, delta = rep(0,p), sigma = par$omega, df = par$nu, log=TRUE  ),silent=TRUE)
	if(is.numeric(v2)) v2=v2 else v2= try(dmvt(r, delta = rep(0,p), sigma = round(par$omega,5), df = par$nu, log=TRUE  ),silent=TRUE)
	if(is.numeric(v2)) v2=v2 else v2= try(dmvt(r, delta = rep(0,p), sigma = round(par$omega,4), df = par$nu, log=TRUE  ),silent=TRUE)
	if(is.numeric(v2)) v2=v2 else v2= try(dmvt(r, delta = rep(0,p), sigma = round(par$omega,3), df = par$nu, log=TRUE  ),silent=TRUE)
	if(is.numeric(v2)) v2=v2 else v2= try(dmvt(r, delta = rep(0,p), sigma = round(par$omega,2), df = par$nu, log=TRUE  ),silent=TRUE)
	if(is.numeric(v2)) v2=v2 else v2= try(dmvt(r, delta = rep(0,p), sigma = round(par$omega,2), df = par$nu, log=TRUE  ),silent=TRUE)
		if(is.numeric(v2)) v2=v2 else v2= try(dmvt(r, delta = rep(0,p), sigma = round(par$omega,1), df = par$nu, log=TRUE  ),silent=TRUE)

	if(is.numeric(v2)) v2=v2 else v2= try(dmvt(r, delta = rep(0,p), sigma = round(par$omega,0), df = par$nu, log=TRUE  ),silent=TRUE)

	v3 = numeric(n)
	for (k in 1:n) v3[k] = log(pmvt2(lower= rep(-Inf,p), upper=q.star[k,], delta = rep(0,p), sigma = round(par$delta,3), df=p+par$nu ))
	
	lval = v1 + v2 + v3
	if (!logd) val = exp(lval)
	else val = lval
	return(val)
}


pmvt2 <- function(lower=NULL, upper= NULL, delta= NULL, sigma= NULL, df=NULL) {
	nu = df
	if (length(lower) == 1) val = pt( (upper-delta)/sqrt(sigma), df=round(nu))
	else val = pmvt(lower=lower, upper=upper, delta=delta, sigma=sigma, df=round(nu))
	return(val)
}

# library(alr3)
# data(ais)
# x=ais[,c("Wt","Bfat","BMI","SSF","Ht")]
# temp=EM(x,2,1)

# Rprof("boot.out")
# temp=EM(x,2,1)
# Rprof(NULL)


# vec1=mvrnorm(100,mu=c(0,0,0,0,0,0,0),Sigma=diag(7))
# vec2=mvrnorm(100,mu=c(20,20,20,20,20,20,20),Sigma=diag(7))
# x=rbind(vec1,vec2)

# test=EM(x,2,1)



# vec1=mvrnorm(100,mu=c(0,0,0,0,0,0,0),Sigma=diag(7))
# vec2=mvrnorm(100,mu=c(20,20,20,20,20,20,20),Sigma=diag(7))
# x=rbind(vec1,vec2)

# test=EM(x,2,1)






################################
################################
##Three functions in here:
#### 1) test.gamma for finding the ML estimate for the grand lambda
#### 2) make.Z for converting a matrix to the Z-scale
#### 3) update.UDV for updating the ideal point estimates
#### *) Small additonal ones at the end 
################################
################################

################################
################################
##1) Finding gamma
#test.gamma.pois_gibbs<-function(gamma.try,Theta.last.0=Theta.last[row.type=="count",],votes.mat.0=votes.mat[row.type=="count",],emp.cdf.0,cutoff.seq=NULL){
test.gamma.pois_gibbs<-function(gamma.try,Theta.last.0,votes.mat.0,emp.cdf.0,cutoff.seq=NULL){
	gamma.try[gamma.try>2]<- 2
	gamma.try[gamma.try< -2]<- -2
	gamma.try<-exp(gamma.try)
	votes.mat<-votes.mat.0
	taus.try<-NULL
	count.seq<-cutoff.seq
	if(length(cutoff.seq)==0) count.seq<-seq(-1,max(votes.mat)+2,1)
	taus.try<-count.seq*0

	analytic.cdf<-count.seq
	a.int<-qnorm(mean(votes.mat==0))
	taus.try[count.seq>=0]<-(a.int+gamma.try[1]*count.seq[count.seq>=0]^(gamma.try[2]))
	taus.try[count.seq<0]<- -Inf
	taus.try<-sort(taus.try)
	taus.try[1]<--Inf	
	find.pnorm<-function(x){
		a<-x[1]
		b<-x[2]
		coefs<-c( -1.82517672, 0.51283415, -0.81377290, -0.02699400, -0.49642787, -0.33379312, -0.24176661, 0.03776971)
		x<-c(1, a, b, a^2,b^2, log(abs(a-b)),log(abs(a-b))^2,a*b)
		(sum(x*coefs))
	}
	lik<-((pnorm(taus.try[votes.mat+2 ]-Theta.last.0)-pnorm(taus.try[votes.mat+1 ]-Theta.last.0)) )
	log.lik<-log(lik)
	
	which.zero<-which(lik==0)
	a0<-taus.try[votes.mat[which.zero]+2]-Theta.last.0[which.zero]
	b0<-taus.try[votes.mat[which.zero]+1]-Theta.last.0[which.zero]
	log.lik[which.zero]<-apply(cbind(a0,b0),1,find.pnorm)
	log.lik[is.infinite(lik)&lik>0]<-max(log.lik[is.finite(lik)])
	log.lik[is.infinite(lik)&lik<0]<-min(log.lik[is.finite(lik)])

	thresh<-min(-1e30, min(log.lik[is.finite(log.lik)],na.rm=TRUE))
	log.lik[log.lik<thresh]<-thresh
	dev.out<--2*sum(log.lik,na.rm=TRUE)+sum(log(gamma.try)^2)
	
	
	return(list("deviance"=dev.out,"tau"=taus.try))

}	



################################
################################
##2) Converting a matrix to the Z scale
	make.Z_gibbs<-function(
			Theta.last.0=Theta.last, 
			votes.mat.0=votes.mat,row.type.0=row.type, 
			n0=n, k0=k, params=NULL, iter.curr=0,empir=NULL,cutoff.seq.0=NULL,missing.mat.0=NULL,lambda.lasso,proposal.sd,scale.sd, max.optim,step.size,maxdim.0,tau.ord.0
			){
		
		tau.ord<-tau.ord.0
		Theta.last<-Theta.last.0;votes.mat<-votes.mat.0;
		row.type<-row.type.0; n<-n0; k<-k0
		cutoff.seq<-cutoff.seq.0
		sigma<-1
		row.type<-row.type.0; votes.mat<-votes.mat.0
		Z.next<-matrix(NA,nrow=n,ncol=k)
		missing.mat<-missing.mat.0
		i.gibbs<-iter.curr
		maxdim<-maxdim.0
		#print(i.gibbs)
		#print(iter.curr)
	

	
			if(sum(row.type=="bin")>2){
				
			Z.next[row.type=="bin",][votes.mat[row.type=="bin",]==1]<-rtruncnorm(sum(votes.mat[row.type=="bin",]==1), mean=sigma^.5*Theta.last[row.type=="bin",][votes.mat[row.type=="bin",]==1], a=0, sd=sigma^.5)
		Z.next[row.type=="bin",][votes.mat[row.type=="bin",]==0]<-rtruncnorm(sum(votes.mat[row.type=="bin",]==0), mean=sigma^.5*Theta.last[row.type=="bin",][votes.mat[row.type=="bin",]==0], b=0 , sd=sigma^.5)
		Z.next[row.type=="bin",][missing.mat[row.type=="bin",]==1]<-rnorm(sum(missing.mat[row.type=="bin",]==1), mean=sigma^.5*Theta.last[row.type=="bin",][missing.mat[row.type=="bin",]==1], sd=sigma^.5)
		
		
		pars.max<-1
		accept.out<-prob.accept<-1

				}

	

		if(sum(row.type=="count")>2){
		#begin type = "pois"
		num.sample.count<-sum(row.type=="count")*k
		accept.out<-prob.accept<-NA

		dev.est<-function(x) test.gamma.pois_gibbs(x,votes.mat.0=votes.mat[row.type=="count",], Theta.last.0=Theta.last[row.type=="count",],cutoff.seq=cutoff.seq)$de
		dev.est1<-function(x) test.gamma.pois_gibbs(c(x,params[2]),votes.mat.0=votes.mat[row.type=="count",], Theta.last.0=Theta.last[row.type=="count",],cutoff.seq=cutoff.seq)$de
		dev.est2<-function(x) test.gamma.pois_gibbs(c(params[1],x),votes.mat.0=votes.mat[row.type=="count",], Theta.last.0=Theta.last[row.type=="count",],cutoff.seq=cutoff.seq)$de
	
	
	#range.opt<-c(params-1,params+1)
	grad.est<-function(x,del=.001){
		grad1<-(dev.est1(x[1]+del)-dev.est1(x[1]-del))/(4*del)
		grad2<-(dev.est2(x[2]+del)-dev.est2(x[2]-del))/(4*del)
		grad.all<-c(grad1,grad2)
		grad.all[!is.finite(grad.all)]<-0
		grad.all
		}
	grad.est.1<-function(x,del=.001){
		#grad1<-(dev.est1(x[1]+del)-dev.est1(x[1]-del))/(4*del)
		#grad2<-0#(dev.est2(x[2]+del)-dev.est2(x[2]-del))/(4*del)
		#c(grad1,grad2)
		grad1<-dev.est1(x[1]+del)-2*dev.est(x)+dev.est1(x[1]-del)
		grad1<-grad1/(2*del)
		grad1[!is.finite(grad1)]<-0
		grad1
		}
	grad.est.2<-function(x,del=.001){
		#grad1<-0#(dev.est1(x[1]+del)-dev.est1(x[1]-del))/(4*del)
		#grad2<-(dev.est2(x[2]+del)-dev.est2(x[2]-del))/(4*del)
		#grad2[!is.finite(grad2)]<-0
		grad2<-dev.est2(x[2]+del)-2*dev.est(x)+dev.est2(x[2]-del)
		grad2<-grad2/(2*del)
		#c(grad1,grad2)
		grad2
		}
	grad.est.12<-function(x,del=.001){
		#grad1<-0#(dev.est1(x[1]+del)-dev.est1(x[1]-del))/(4*del)
		#grad12<-(dev.est(x+del)-dev.est(x-del))/(8*del)
		grad12<-dev.est(x+del)-dev.est(x+del*c(1,-1))-dev.est(x+del*c(-1,1))+dev.est(x-del)
		grad12<-grad12/(4*del)
		grad12[!is.finite(grad12)]<-0
		grad12
		#c(grad1,grad2)
		}


			##Hamiltonian step
		##Code adapted from Neal, Handbook of MCMC
		gamma.opt<-"Hamiltonian Step"
		step.size<-min(step.size,2)
		
		#print("step size")
		#print(step.size)
		
		q0<-q<-params
		p0<-p<-rnorm(length(q),0,1)
		current_U<-dev.est(q0)/2

		num.HMC<-10
		for(i.HMC in 1:(num.HMC)){
		if(i.HMC%%3==1){
		hess.mat<-matrix(NA,nrow=2,ncol=2)
		hess.mat[1,1]<-grad.est.1(q)
		hess.mat[1,2]<-hess.mat[2,1]<-grad.est.12(q)
		hess.mat[2,2]<-grad.est.2(q)
		hess.mat[!is.finite(hess.mat)]<-0
		hess.mat<-.5*(hess.mat+t(hess.mat))
		epsilon<-ginv(hess.mat)*step.size/10
		}
		q<-q+epsilon%*%p
		p<-p-epsilon%*%grad.est(q)
		if(abs(dev.est(q)/2 - current_U)>10000) {revert.q<-TRUE;q<-q0;break}
		revert.q<-FALSE
		}
		q<-q+epsilon%*%p
		p<-p-epsilon%*%grad.est(q)/2
		p<- -p
		
		if(revert.q) q<-q0
		
		current_K<-sum(p0^2)/2
		proposed_U<-dev.est(q)/2
		proposed_K<-sum(p^2)/2
		proposal.dev<-proposed_U+proposed_K
		current.dev<-current_U+current_K
		gamma.cand<-q
		#print("Log Probs")
		#print(c(proposal.dev,current.dev))
		if(exp(current.dev-proposal.dev)>runif(1)){
			##Do if accept
			pars.max<-gamma.next<-(gamma.cand)
			#print("accept")
			accept.out<-prob.accept<-1
			} else {
			##Do if reject
			pars.max<-gamma.next<-(params)
			#print("reject")
			accept.out<-prob.accept<-0
			}
		
		#print("Proposal sd")
		#print(proposal.sd)
		taus<-test.gamma.pois_gibbs(gamma.next,votes.mat.0=votes.mat[row.type=="count",], Theta.last.0=Theta.last[row.type=="count",],cutoff.seq=cutoff.seq)$tau
		taus<-sort(taus)
		taus[1]<--Inf
		#print("Range of tau")
		#print(range.opt)
		#print(gamma.opt)
		#print(range(taus))
		#print("Beta parameters")
		#print(pars.max)
		
		Z.next.temp<-rtruncnorm(num.sample.count, mean=Theta.last[row.type=="count",], a=taus[votes.mat[row.type=="count",]+1 ],b=taus[votes.mat[row.type=="count",]+2 ] )

		
		Z.next[row.type=="count",]<-matrix(Z.next.temp,nrow=sum(row.type=="count"))



		}
		
	return(list("Z.next"=Z.next,"params"=(pars.max),"accept"=accept.out,"prob"=prob.accept,"proposal.sd"=proposal.sd,"step.size"=step.size,"tau.ord"=tau.ord))
	}##Closes out make.Z function


################################
################################
##2) Converting a matrix to the Z scale
	update_UDV_gibbs<-function(
	Z.next.0=Z.next,
	k0=k, n0=n, lambda.lasso.0=lambda.lasso,lambda.shrink.0=lambda.shrink,
	Dtau.0=Dtau,
	votes.mat.0=votes.mat, iter.curr=0,row.type.0,missing.mat.0=missing.mat,maxdim.0,
	V.last
	){

	missing.mat<-missing.mat.0
	Dtau<-Dtau.0;
	Z.next<-Z.next.0;
	votes.mat<-votes.mat.0;
	n<-n0; k<-k0; 
	lambda.lasso<-lambda.lasso.0
	row.type <- row.type.0
	maxdim<-maxdim.0

	#Declare some vectors
	sigma<-1
	ones.r<-rep(1,k0)
	ones.c<-rep(1,n0)

	#Update intercepts
	mu.r<-rowMeans(Z.next)#-ones.c%*%t(mu.c)-Theta.last.0+mu.grand)
	mu.r<-mu.r*n/(n+1)+rnorm(length(mu.r),sd=1/k)

	mu.c<-colMeans(Z.next)#-mu.r%*%t(ones.r)-Theta.last.0+mu.grand)
	mu.c<-mu.c*k/(k+1)+rnorm(length(mu.c),sd=1/n)
	mu.grand<-mean(Z.next)
	mu.grand<-mu.grand*(n*k)/(n*k+1)+rnorm(1,sd=1/(n*k))
	
	mean.mat<-ones.c%*%t(mu.c)+mu.r%*%t(ones.r)-mu.grand
	
	if(length(unique(row.type))>1){
		is.bin<-sum(row.type=="bin")>0
		is.count<-sum(row.type=="count")>0
		is.ord<-sum(row.type=="ord")>0
		
		#print("Two Means Being Used")
		mean.c.mat<-matrix(NA,nrow=n,ncol=k)
		if(is.bin) mu.c.bin<-colMeans(Z.next[row.type=="bin",])
		if(is.count) mu.c.count<-colMeans(Z.next[row.type=="count",])
		if(is.ord) mu.c.ord<-colMeans(Z.next[row.type=="ord",])

		if(is.bin) mu.c.bin<-mu.c.bin*k/(k+1)+rnorm(length(mu.c.bin),sd=1/sum(row.type=="bin"))
		if(is.count) mu.c.count<-mu.c.count*k/(k+1)+rnorm(length(mu.c.count),sd=1/sum(row.type=="count"))
		if(is.ord) mu.c.ord<-mu.c.ord*k/(k+1)+rnorm(length(mu.c.ord),sd=1/sum(row.type=="ord"))

		if(is.bin) mean.c.mat[row.type=="bin",]<-ones.c[row.type=="bin"]%*%t(mu.c.bin)
		if(is.count) mean.c.mat[row.type=="count",]<-ones.c[row.type=="count"]%*%t(mu.c.count)
		if(is.ord) mean.c.mat[row.type=="ord",]<-ones.c[row.type=="ord"]%*%t(mu.c.ord)

	
		mean.grand.mat<-matrix(NA,nrow=n,ncol=k)
		if(is.bin) mean.grand.mat[row.type=="bin",]<-mean(Z.next[row.type=="bin",])+rnorm(1,sd=1/(sum(row.type=="bin")*k))
		if(is.count) mean.grand.mat[row.type=="count",]<-mean(Z.next[row.type=="count",])+rnorm(1,sd=1/(sum(row.type=="count")*k))
		if(is.ord) mean.grand.mat[row.type=="ord",]<-mean(Z.next[row.type=="ord",])* (sum(row.type=="count")*k)/(sum(row.type=="ord")*k+1)+rnorm(1,sd=1/(sum(row.type=="ord")*k))


	mean.mat<-mean.c.mat+mu.r%*%t(ones.r)-mean.grand.mat
	
	}

	Z.starstar<-svd.mat<- Z.next- mean.mat



	

	#Take svd, give each column an sd of 1 (rather than norm of 1)
	num.zeroes<-1#colMeans(votes.mat!=0)
	svd.mat.0<-svd.mat
	svd.mat.0<-apply(svd.mat.0,2,FUN=function(x) x-mean(x))
	drop.cols<-colMeans(svd.mat.0^2)^.5<1e-2
	drop.cols[!is.finite(drop.cols)]<-TRUE
	if(sum(drop.cols)>0) svd.mat[,drop.cols]<-rnorm(nrow(svd.mat)*sum(drop.cols))/3
	
	save(svd.mat,file="svd.mat")
	
	svd.mat[is.na(svd.mat)]<-0
	svd.mat.dum<-svd.mat
	
	wts.dum<-rep(1,nrow(svd.mat))
	wts.dum[row.type=="bin"]<-1/sum(1-missing.mat[row.type=="bin",])^.5
	wts.dum[row.type=="count"]<-1/sum(1-missing.mat[row.type=="count",])^.5
	wts.dum[row.type=="ord"]<-1/sum(1-missing.mat[row.type=="ord",])^.5


	wts.dum<-wts.dum/mean(wts.dum)

	svd.dum<-svd(svd.mat*wts.dum,nu=maxdim,nv=maxdim)
	#svd.dum<-irlba(svd.mat*wts.dum,nu=maxdim,nv=maxdim,V=V.last)
	svd.dum$u[!is.finite(svd.dum$u)]<-0
	svd.dum$d[!is.finite(svd.dum$d)]<-0
	svd.dum$v[!is.finite(svd.dum$v)]<-0


	svd.dum$d<-(t(svd.dum$u)%*%svd.mat%*%svd.dum$v)
	which.rows<-which(rowMeans(svd.dum$d^2)^.5<1e-4)
	which.cols<-which(colMeans(svd.dum$d^2)^.5<1e-4)
	svd.dum$d[which.rows,]<-rnorm(length(svd.dum$d[which.rows,]),sd=.001)
	svd.dum$d[which.cols,]<-rnorm(length(svd.dum$d[which.cols,]),sd=.001)
	svd2<-svd(svd.dum$d,nu=maxdim,nv=maxdim)
	#print(dim(svd.dum$d))
	#svd2<-irlba(svd.dum$d,nu=maxdim,nv=maxdim)
	svd2$u[is.na(svd2$u)|is.infinite(svd2$u)]<-0
	svd2$d[is.na(svd2$d)|is.infinite(svd2$d)]<-0
	svd2$v[is.na(svd2$v)|is.infinite(svd2$v)]<-0
	
	svd.dum$u<-svd.dum$u%*%(svd2$u)
	svd.dum$v<-svd.dum$v%*%(svd2$v)

	svd.dum$d<-svd2$d
	
	svd0<-svd.dum
	
	svd0$v<-t(t(svd0$u)%*%svd.mat.0)
	svd0$v<-apply(svd0$v,2,FUN=function(x) x/sum(x^2)^.5)
	svd0$d<-diag(t(svd0$u)%*%svd.mat.0%*%svd0$v)
	sort.ord<-sort(svd0$d,ind=T,decreasing=T)$ix
	svd0$u<-svd0$u[,sort.ord]
	svd0$v<-svd0$v[,sort.ord]
	svd0$d<-svd0$d[sort.ord]
	svd0$u<-svd0$u*(n-1)^.5
	svd0$v<-svd0$v*(k-1)^.5
	svd0$d<-svd0$d*((n-1)*(k-1))^-.5
	Theta.last.0<-svd.mat
	Theta.last<-Theta.last.0+(ones.c%*%t(mu.c)+mu.r%*%t(ones.r)-mu.grand)

	#Update d; follows from Blasso and DvD
	Y.tilde<-as.vector(svd.mat)
	if(n>length(Dtau)) Dtau[(length(Dtau)+1):n]<-1
	A<- (n*k)*diag(n)+diag(as.vector(Dtau^(-1)))
	#if(n>k) for(i.A in (k+1):n) A[i.A,i.A]<-(Dtau^-1)[i.A]*0+1e20
	#gA<-ginv(A)
	gA<-A*0+NA
	gA[1:maxdim,1:maxdim]<-ginv(A[1:maxdim,1:maxdim])
	gA[is.na(gA)]<-0
	XprimeY<-sapply(1:maxdim, FUN=function(i, svd2=svd0, Z.use=svd.mat) sum(Z.use*(svd2$u[,i]%*%t(svd2$v[,i]))))
	if(length(XprimeY)<dim(gA)[1]) XprimeY[(length(XprimeY)+1) :dim(gA)[1]]<-0
	D.post.mean<- as.vector(gA%*%XprimeY)
	D.post.var.2<-gA

##Sample D and reconstruct theta, putting intercepts back in
	D.post<-rep(NA,length(D.post.mean))
	D.post[1:maxdim]<-as.vector(mvrnorm(1, mu=D.post.mean[1:maxdim], D.post.var.2[1:maxdim,1:maxdim] ) )/sigma^.5
	D.post[is.na(D.post)]<-0
	abs.D.post<-abs(D.post)

##Calculate MAP and mean estimate
	D.trunc<-pmax(svd0$d-lambda.lasso.0,0)
	
	U.last<-svd0$u
	V.last<-svd0$v

##Update U, V
	prior.var<-.25	
	U.last<-t(t(U.last)*(D.post[1:maxdim]^2/(D.post[1:maxdim]^2+1/(4*k))))
	ran.mat<-sapply(1/(D.post[1:maxdim]^2*(n-1)+.25),FUN=function(x) rnorm(nrow(U.last),mean=0,sd=x^.5))
	U.last<-U.last+ran.mat
	V.last<-t(t(V.last)*D.post[1:maxdim]^2/(D.post[1:maxdim]^2+1/(4*n)))
	ran.mat<-sapply(1/(D.post[1:maxdim]^2*(k-1)+.25 ),FUN=function(x) rnorm(nrow(V.last),mean=0,sd=x^.5))
	V.last<-V.last+ran.mat

	U.next<-U.last
	V.next<-V.last
	
	#Uncommenting the next two lines eliminates the U and V
	#Gibbs update.
	#U.next<-svd0$u
	#V.next<-svd0$v

	Theta.last.0<-U.next%*%diag(D.post[1:maxdim])%*%t(V.next)
	Theta.last<-Theta.last.0+mean.mat
	Theta.last[row.type=="bin",]<-Theta.last[row.type=="bin"]-mean(Theta.last[row.type=="bin",][missing.mat[row.type=="bin"]==0])+ qnorm(mean(votes.mat[row.type=="bin",][missing.mat[row.type=="bin"]==0]))
	#Not sure what the next three lines did.
	#chisq.ran<-rchisq(1,n*k)#+rchisq(maxdim,1)
	#sigma<-(sum((Z.next-Theta.last)^2))/sum(chisq.ran)
	#Theta.last<-Theta.last
	Theta.mode<- U.next%*%diag(D.trunc[1:maxdim])%*%t(V.next)+mean.mat
	##Update muprime, invTau2, lambda.lasso
	muprime<-(abs(lambda.lasso*sqrt(sigma)/(abs.D.post)))
	invTau2<-sapply(1:maxdim, FUN=function(i) rinv.gaussian(1, muprime[i], (lambda.lasso^2 ) ) )
	Dtau<-abs(1/invTau2)
	lambda.lasso<-rgamma(1, shape=maxdim+1  , rate=sum(Dtau[1:maxdim])/2+1.78  )^.5
	lambda.shrink<-(maxdim/(sum(Dtau[1:maxdim])/2+1.78))^.5
	
	
	return(list(
	"Theta.last"=Theta.last,
	"U.next"=U.next,
	"V.next"=V.next,
	"lambda.lasso"=lambda.lasso,
	"lambda.shrink"=lambda.shrink,
	"D.trunc"=D.trunc,
	"D.post"=D.post,
	"Theta.mode"=Theta.mode,
	"svd0"=svd0, 
	"Dtau"=Dtau,
	"D.ols"=svd0$d
	))
}



	expit<-function(x) exp(x)/(1+exp(x))

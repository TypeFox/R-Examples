
########################################
########################################
####The functions
########################################
########################################

#Takes type linear, tobit, probit

#scale.type: none, expand.treat, expand.X

sparsereg<-function(y, X, treat=NULL, EM=FALSE, gibbs=200, burnin=200, thin=10,  type="linear",scale.type="none", baseline.vec=NULL, id=NULL, id2=NULL, id3=NULL, save.temp=FALSE, conservative=TRUE){

EM.all<-FALSE
if(EM) EM.all<-TRUE
#To be added in later!
#type="linear"
#trunc=NULL 
lower.trunc=FALSE
one.constraint=FALSE
asymp<-TRUE
if(EM){
	thin=2
	burnin=0
}
	
	##Warnings
	if(!scale.type%in%c("none","TX","TTX","TT","XX")) stop("scale.type should be one of none, TX, TTX, TT, or XX")
	if(!type%in%c("linear","probit")) stop("type should be one of linear or probit")
	#if(type==("probit")&EM==TRUE) stop("EM is not currently implemented for a sparse probit regression.\n Please try the Gibbs implementation.")
	if(sum(!is.finite(y))>0) stop("Please remove missing or infinite values from y")
	if(sum(!is.finite(X))>0) stop("Please remove missing or infinite values from X")
	#if(nrow(as.matrix(X))!=length(y)) stop("y and X need to have the same length")
	if(length(id)>0& length(id)!=length(y)) stop("id and y need to have the same length")
	if(length(id2)>0& length(id2)!=length(y)) stop("id2 and y need to have the same length")
	if(length(id3)>0& length(id3)!=length(y)) stop("id3 and y need to have the same length")

	###################################
	###################################
	#### Declare variables
	###################################
	###################################
	#set.group<-length(group)!=0
	
	#if(length(X)>0) if(length(group)==0) group<-ifelse(length(X)>0,rep(1,ncol(X)),999)
	#if(length(group)!=0) group<-as.numeric(as.factor(group))
	
	id<-as.numeric(as.factor(id))
	id2<-as.numeric(as.factor(id2))
	id3<-as.numeric(as.factor(id3))

	if(length(treat)!=0) treat<-data.frame(treat)
	n<-length(y)

if(!EM) cat("Step 1 of 4: Formatting Data","\n")
if(EM)  cat("Step 1 of 3: Formatting Data","\n")
	if(length(X)>0){
		if(length(colnames(X))==0) colnames(X)<-paste("X",1:ncol(X),sep="")
	X<-X[,apply(X,2,sd)>0]
	}
	
	if(type=="tobit") y0<-y
	
	if(type=="probit") {
	y0<-(y>mean(y))*1
	y[y0==1]<-dnorm(0)/pnorm(0)
	y[y0==0]<--dnorm(0)/pnorm(0)
	
	
	probit.int<-0#mean(y)
	y<-y-probit.int
	}
	
	y<-as.vector(y)
	if(length(X)>0){
		X<-as.matrix(X)
		X<-apply(X,2,as.numeric)
		}
	if(length(treat)==0) X.big<-X.lin<-X
	drop.string<-"IPREFERONEDIRECTIONTOLEDZEPPELIN"
	if(length(treat)!=0){
		treat<-as.matrix(treat)
		if(length(baseline.vec)>0){
		for(i.b in 1:length(baseline.vec)) 
		treat[treat[,i.b]==baseline.vec[i.b],i.b]<-drop.string
		}
			
	dummy.col<-function(t1){
		t1<-as.matrix(t1)
		if(length(colnames(t1))==0) colnames(t1)<-"treat"
		out<-sapply(sort(unique(t1)),FUN=function(x) x==t1)*1
		colnames(out)<-paste(colnames(t1)[1],sort(unique(t1)), sep="_")
		invisible(out)
	}
	if(length(treat)>0) {
	treat0<-as.matrix(treat)
	treat.mat.0<-NULL
	for(i.t in 1:ncol(treat)) treat.mat.0<-cbind(treat.mat.0,dummy.col(treat0[,i.t]))
	}
	if(length(X)>0) X.lin<-cbind(X,treat.mat.0) else X.lin<-treat.mat.0

	X.lin<-apply(X.lin,2,as.numeric)
	colnames(X.lin)<-gsub(":","_",colnames(X.lin))
	}
	drop.cols<-grep(drop.string,colnames(X.lin))
	if(length(drop.cols)>0) X.lin<-X.lin[,-drop.cols]
	
	
if(scale.type=="none"){
	if(length(X)>0) c.cluster<-c(rep(1,ncol(X)),rep(2,ncol(X.lin)-ncol(X))) else
		c.cluster<-rep(1,ncol(X.lin))
	X.lin<-apply(X.lin,2,FUN=function(x) x-mean(x) )
	X.big<-X.lin
	if(one.constraint) c.cluster<-rep(1,ncol(X.big))
	}


	if(scale.type=="TTX"){
			treat<-as.matrix(treat)
			all.vars<-make.threewayinter(X,treat,treat)
			X.big<- all.vars$big.X
			c.cluster<-all.vars$c.clust
	}


	if(scale.type=="TX"){
			if(length(X)==0) X<-matrix(1,nrow=length(y),ncol=1)
			treat<-as.matrix(treat)
			all.vars<-make.inter(X,treat)
			X.big<- all.vars$X
			c.cluster<-all.vars$c.clust

	}

	if(scale.type=="TT"&(length(X)==0)){
			treat<-as.matrix(treat)
			X.lin<-apply(X.lin,2,as.numeric)
			all.vars<-make.inter(X.lin,X.lin)
			X.big<- all.vars$X
			c.cluster<-all.vars$c.clust

	}


	if(scale.type=="TT"&(length(X)>0)){
			all.vars<-make.inter.treat(X,treat)
			X.big<- all.vars$X
			c.cluster<-all.vars$c.clust
		}

	if(scale.type=="XX"){
			all.vars<-make.interXX(X,treat)
			X.big<- all.vars$X
			c.cluster<-all.vars$c.clust
		}


	drop.cols<-grep(drop.string,colnames(X.big))
	if(length(drop.cols)>0) X.big<-X.big[,-drop.cols]
	
	###################################
	###################################
	#### Eliminate correlated vars from X
	###################################
	###################################
	X.big[abs(X.big)<1e-10]<-0

	X.big<-apply(X.big,2,FUN=function(x) x-mean(x))
	keeps<-apply(X.big,2,sd)>0&(apply(X.big,2,FUN=function(x) sum(is.na(x)))==0)

	#if(set.group)	{c.cluster<-group
	#if(length(group)!=ncol(X.big)) stop("group should have one entry for every element of the full X matrix")
	#}

	X.big<-X.big[,keeps]
	c.cluster<-c.cluster[keeps]

	
	cor1<-cor(X.big)
	diag(cor1)<-0
	cors.run<-NULL
for(i in 1:(ncol(X.big)-1)) for(j in (i+1):ncol(X.big)){
	if(cor1[i,j]^2>0.999) cors.run<-rbind(cors.run,c(i,j))
}

	if(length(cors.run)>0){
	drops<-as.vector(cors.run[,2])
	X.big<-X.big[,-drops]
	c.cluster<-c.cluster[-drops]
	}
	

	scale.big<-sapply(colnames(X.lin),grep,colnames(X.big))
	matrix.convert<-matrix(0,nrow=length(colnames(X.lin)),ncol=ncol(X.big))
	rownames(matrix.convert)<-colnames(X.lin)
	colnames(matrix.convert)<-colnames(X.big)
	for(i.scale in 1:ncol(X.lin)) matrix.convert[i.scale,unlist(scale.big[i.scale])]<-1

	#if(length(group)>0){if(group=="bylevel") c.cluster<-as.numeric(as.factor(colSums(matrix.convert)))}

	scale.back<-apply(X.big,2,FUN=function(x) sd(x))

	for(i.scale in 1:ncol(X.big)) X.big[,i.scale]<-X.big[,i.scale]/scale.back[i.scale]
	
	scale.y<-sd(y)
	y<-(y-mean(y))/scale.y

	n.cluster<-length(unique(c.cluster))


	###################################
	###################################
	#### Format containers; initialize
	###################################
	###################################

	n<-length(y)
	p<-ncol(X.big)
	c.cluster<-c.cluster[1:p]

	beta.mean<-matrix(NA,nrow=p,ncol=gibbs)
	beta.mode<-matrix(NA,nrow=p,ncol=gibbs)
	lambda.run<-matrix(NA,nrow=p,ncol=gibbs)
	beta.ci<-matrix(NA,nrow=p,ncol=gibbs)

	sigma.sq.run<-rep(NA,gibbs)
	rownames(beta.mean)<-rownames(beta.mode)<-colnames(X.big)
	beta.curr.mode<-beta.curr<-rnorm(p,sd=.1)
	

	ps.sigmasq<-sigma.sq<-mean((y-X.big%*%beta.curr)^2)
	if(type=="probit") sigma.sq<-1

	##For Dirichlet Dirichlet Laplacian; not used anymore
	psi<-phi<-rep(1/p,p)
	grand.lambda<-grandTau<-1
	lambda.shrink<-lambda.all<-rep(0,p)
	rel.wts<-1
	alpha0<-.25
	
	for(i in unique(c.cluster)) {
		lambda.all[c.cluster==i]<-(2*sum(c.cluster==i)/sum(abs(beta.curr)[c.cluster==i]))^.5
		}
	lambda.shrink<-lambda.all
	if(one.constraint) lambda.all[1:p]<-(2*p/sum(abs(beta.curr)))^.5
	k.vec<-rep(NA,length(c.cluster))
	for(i in 1:n.cluster) k.vec[c.cluster==i]<-sum(c.cluster==i)
	#if(ncol(X.big)<150) block<-TRUE
	#block<-FALSE
	#if(block) XprimeX<-t(X.big)%*%X.big


	##Initialize Random Effects
	ran.effect<-rep(0,n)
	sigma.sq.b<-1
	k.raneff<- length(unique(id))+length(unique(id2))+length(unique(id3))

	if(k.raneff>0){
	ran.effect1<-ran.effect2<-ran.effect3<-rep(0,n)
	sigma.sq.b<-sigma.sq.b2<-sigma.sq.b3<-1

	re1<-updateREs(y, fits.main=as.vector(X.big%*%beta.curr), ran.effect1, ran.effect2, ran.effect3, id, id2, id3,sigma.sq.b, sigma.sq.b2,sigma.sq.b3,sigma.sq=1,fix.eff=FALSE,EM0=EM)
	sigma.sq.b<-re1$sigmas[1]
	sigma.sq.b2<-re1$sigmas[2]
	sigma.sq.b3<-re1$sigmas[3]
	ran.effect1<-re1$REs[,1]
	ran.effect2<-re1$REs[,2]
	ran.effect3<-re1$REs[,3]
	ran.effect<-rowSums(re1$REs)
	} 

if(!EM) cat("Step 2 of 4: Beginning Burnin Period","\n")
if(EM) cat("Step 2 of 3: Beginning EM","\n")

	beta.conv<-beta.curr
	##Begin gibbs loop
	for(i.gibbs in 1:(burnin+gibbs)){	
			if(i.gibbs>1){
		if(EM.all) {
			if(max(abs(beta.conv-beta.curr))<1e-4) break 
			} else{
			EM<-FALSE
				}
		beta.conv<-beta.curr
	}	
	for(i.thin in 1:thin){#Start thin loop
	muprime.all<-abs(lambda.all*sqrt(sigma.sq))/abs(beta.curr)
	muprime.all[is.na(muprime.all)]<-1e6
	lambda.prime<-lambda.all^2
	if(!EM) {
		invTau2<-sapply(1:length(muprime.all), FUN=function(i) rinv.gaussian(1, muprime.all[i], (lambda.prime[i] ) ) )
	Dtau<-abs(1/invTau2)
	} else{
		Dtau<-NULL
		tau.try<-seq(.001,500,.01)
		for(i.tau in 1:length(muprime.all)) {
			range.tau<-NULL
			#range.tau<-range(rinv.gaussian(100,muprime.all[i.tau],lambda.prime[i.tau]))
			range.tau[1]<-muprime.all[i.tau]/5#,lambda.prime[i.tau])
			range.tau[2]<-muprime.all[i.tau]+muprime.all[i.tau]^1.5/lambda.prime[i.tau]^.5*5
			if(sum(is.na(range.tau))>0) range.tau<-c(100,2000)
			tau.try<-seq(range.tau[1],range.tau[2],length=500)
			probs.tau<-dinv.gaussian(tau.try,muprime.all[i.tau],lambda.prime[i.tau])	
			probs.tau<-probs.tau/sum(probs.tau,na.rm=TRUE)
			Dtau[i.tau]<-sum(probs.tau/tau.try,na.rm=TRUE)
			}
		
		}
		
		


	up1<-updatebeta_cpp(X0=X.big, y0=as.matrix(y-ran.effect), betacurr0=as.matrix(beta.curr), betamode0=as.matrix(beta.curr.mode), lambdavec0=as.matrix(lambda.all*rel.wts), dtau0=as.matrix(Dtau*rel.wts^2), sigmasq0=as.matrix(sigma.sq),ps_sigmasq0=as.matrix(ps.sigmasq),lambdashrink0=as.matrix(lambda.shrink*rel.wts),
	k0=as.matrix(1)
	)

	beta.curr<-up1$beta.mean
	beta.curr.mode<-up1$beta.mode
	beta.curr.ci<-up1$beta.ci
	if(EM) {
		beta.curr<-up1$beta.EM*.5+beta.curr*.5
		beta.curr.ci<-beta.curr
	}
	beta.curr[abs(beta.curr)<1e-10]<-1e-10


##Calculate approximate confidence interval
#if(i.thin==thin){
if(i.thin==thin&!EM){
#	ps.sigmasq<-rinvgamma(1,shape=ps.sigma.sq.shape,scale=sigma.sq.scale)

	for(i.var in 1:p){
        beta.ls<-up1$beta.ols[i.var]
        lambda.use<-lambda.all[i.var]
	in.phi<-abs(n^.5*abs(beta.ls/sigma.sq^.5)-lambda.shrink[i.var]*ps.sigmasq^.5/(sigma.sq^.5*(n-1)) ) 
	p.z<-pnorm(in.phi)
        beta.2<-beta.curr.mode[i.var]^2
	var.beta2<-sigma.sq/(n-1)
        var.betals<-beta.curr.ci[i.var]
	var.lambda<-1/(4*grand.lambda^2)*(p*n^.5+1)/(sum(Dtau*rel.wts^2)/2+1.78)^2+grand.lambda^2*wts.lam[i.var]^4*var.u[i.var]  #sum(Dtau*rel.wts^2)/2+1.78)^.5
        beta.sigma<-sigma.sq.scale#(beta.curr.ci[i.var]*(n-1)+sum(beta.curr^2/Dtau))/2
        alpha.sigma<-ps.sigma.sq.shape#(n^.5-1)/2
        #Variance of pseudo sigma
        var.sigma<-1/(4*ps.sigmasq)*beta.sigma/((alpha.sigma-1)^2*(alpha.sigma-2))
        var.sigma.eps<-1/(4*sigma.sq)*sigma.sq.scale/((sigma.sq.shape-1)^2*(sigma.sq.shape-2))

	

        beta.curr.ci[i.var]<-
                (
                p.z^2*var.beta2+
                n/sigma.sq*(beta.2*dnorm(in.phi)*1/sigma.sq^.5)^2*var.betals+
                n/sigma.sq*(beta.2*dnorm(in.phi)*beta.ls*lambda.shrink[i.var]/(n-1))^2*var.sigma+
                n/var.betals*(beta.2*dnorm(in.phi)*ps.sigmasq/(n-1))^2*var.lambda+
                n*(beta.2*dnorm(in.phi)*((abs(beta.2)/sigma.sq)-lambda.shrink[i.var]*ps.sigmasq/sigma.sq/n))^2*var.sigma.eps
                
                )#*(1.43)^2#qcauchy(.1)/qnorm(.1)=3.84#qt(.05,3)/1.6444=1.43 bc largest t w finite variance
		beta.curr.ci[i.var]<-max(beta.curr.ci[i.var],up1$var.beta[i.var])
	
	}
	
	y.temp<-y-ran.effect

	if(!exists("fits")) fits<-X.big%*%up1$beta.curr
	si<-(y.temp-ran.effect-fits)/n
	ki<-1/n
	df.satter<-sum(si^2)^2/sum(si^4)

	beta.curr.ci<-abs(beta.curr.ci)^.5/rgamma(1,df.satter/2,df.satter/2)^.5

	mean.temp<-up1$beta.ols#beta.curr.mode
	set.temp<-beta.curr.mode==0
	if(sum(set.temp)>0)
        beta.curr.ci[set.temp]<-rnorm(sum(set.temp>0),mean.temp[set.temp],sd=beta.curr.ci[set.temp])
	set.temp<-beta.curr.mode>0
	if(sum(set.temp)>0)
	beta.curr.ci[set.temp]<-rtnorm(sum(set.temp),mean=mean.temp[set.temp],sd=beta.curr.ci[set.temp],
	lower=0,upper=Inf)
	set.temp<-beta.curr.mode<0
	if(sum(set.temp)>0)
	beta.curr.ci[set.temp]<-rtnorm(sum(set.temp),mean=mean.temp[set.temp],sd=beta.curr.ci[set.temp],
	lower=-Inf,upper=0)
	}

	
	if(k.raneff>0){
	re1<-updateREs(y, fits.main=as.vector(X.big%*%beta.curr), ran.effect1, ran.effect2, ran.effect3, id, id2, id3,	sigma.sq.b, sigma.sq.b2, sigma.sq.b3, 1,fix.eff=FALSE,EM0=EM)
	sigma.sq.b<-re1$sigmas[1]
	sigma.sq.b2<-re1$sigmas[2]
	sigma.sq.b3<-re1$sigmas[3]
	ran.effect1<-re1$REs[,1]
	ran.effect2<-re1$REs[,2]
	ran.effect3<-re1$REs[,3]
	ran.effect<-rowSums(re1$REs)
	} 
	
	fits<-as.vector(X.big%*%beta.curr)+ran.effect
	##Update sigmas, lambda
	sigma.sq.shape<-(n-1)/2+p/2+k.raneff/2
	ps.sigma.sq.shape<-(n^.5-1)/2+p/2+k.raneff/2

	if(conservative==FALSE)
	ps.sigma.sq.shape<-(n^(.25)-1)/2+p/2+k.raneff/2
	sigma.sq.scale<-sum((y-fits)^2)/2+ sum(beta.curr^2/(Dtau*rel.wts^2))/2

	sigma.sq<-rinvgamma(1,shape=sigma.sq.shape,scale=sigma.sq.scale)

	ps.sigmasq<-rinvgamma(1,shape=ps.sigma.sq.shape,scale=sigma.sq.scale)

	
	if(EM){
		sigma.sq<-sigma.sq.scale/(sigma.sq.shape-1)
		ps.sigmasq<-sigma.sq.scale/(ps.sigma.sq.shape-1)
	}
	if(type=="probit") sigma.sq<-1


	p.adj<-n^.5
	##Change
	if(conservative==FALSE)
	p.adj<-n^(.5-1/8)
#	p.adj<-n^.75
	lambda.shrink<-lambda.all


	rel.wts<-1
#	grand.lambda<-rgamma(1, shape=p*p.adj+1, rate=sum(Dtau*rel.wts^2)/2+1.78)^.5
##Adjustment after first draft--!!!

	if(EM) {

		lambda.range<-NULL
		lambda.range[1]<-qgamma(0.001, shape=p*p.adj, rate=sum(Dtau*rel.wts^2)/2+1)
		lambda.range[2]<-qgamma(1-0.001, shape=p*p.adj, rate=sum(Dtau*rel.wts^2)/2+1)
		lambda.try<-seq(lambda.range[1],lambda.range[2],length=500)
		lambda.prob<-dgamma(lambda.try, shape=p*p.adj, rate=sum(Dtau*rel.wts^2)/2+1)
		lambda.prob[is.na(lambda.prob)]<-0
		lambda.prob<-lambda.prob/sum(lambda.prob)
		grand.lambda<-sum(lambda.prob*lambda.try^.5)
		#grand.lambda<-p*p.adj/(sum(Dtau*rel.wts^2)/2+1)
	} else{

		grand.lambda<-rgamma(1, shape=p*p.adj, rate=sum(Dtau*rel.wts^2)/2+1)^.5
}
		wts.lam<-NULL
		#for(i.g in 1:p) wts.lam[i.g]<-rgig(1, .5, 1,abs(beta.curr/sigma.sq^.5)[i.g])
		#wts.lam<-rgamma(p,p.adj+1/2,(abs(beta.curr)/sqrt(sigma.sq))^.25)
		#wts.lam<-wts.lam/mean(wts.lam)
		#print(beta.curr[1:10])	

	#var.u[whatev]<- gets the estimated variance

		u.calc<-function(beta,alpha.func=alpha0,EM0=EM){
			u.try<-seq(.025,20,.025)
			#p.al<-exp(-u.try*beta-u.try^-alpha.func)*u.try^(-(2+n^.25))
			p.al.exp<-(-u.try*beta-u.try^-alpha.func)-(2+n^.25)*(log(u.try))
			p.al.exp<-p.al.exp-max(p.al.exp)
			p.al.exp[p.al.exp< -20]<- -20
			p.al<-exp(p.al.exp)
			
			p.al<-p.al/sum(p.al)
			
			var.out<-sum(u.try^2*p.al)-(sum(u.try*p.al))^2
			if(!EM0) u.samp<-sample(size=1,u.try,prob=p.al)
			if(EM0) u.samp<-sum(p.al*u.try)
			#print(which(out==u.try)/length(u.try))
			out<-list("var.out"=var.out,"u.samp"=u.samp)
			#sum(u.try*p.al)/sum(p.al)
			return(out)
			}
			
		var.u<-NULL
			
		#wts.lam<-sapply(grand.lambda*abs(beta.curr/sigma.sq^.5),u.calc)

		for(i.u in 1:length(beta.curr)) {
			temp.calc<-u.calc(grand.lambda*abs(beta.curr[i.u]/sigma.sq^.5))
			wts.lam[i.u]<-temp.calc$u.samp
			var.u[i.u]<-temp.calc$var.out
		}


		alpha.calc<-function(x){
			-sum(abs(wts.lam)^x)+2*sum(log(wts.lam))+    p*log(x)-p*lgamma(4/x)-x	}
		
	#print("here?")
		alpha.try<-seq(.025,10,.025)
		log.alpha<-sapply(alpha.try,alpha.calc)
		log.alpha<-log.alpha-max(log.alpha)
		log.alpha[log.alpha< -20]<--20
		#print(log.alpha)
		wt.alpha<-exp(log.alpha)#/sum(exp(log.alpha))
	
	
		if(!EM) alpha0<-sample(size=1,alpha.try,prob=wt.alpha)
		if(EM) {
			wt.alpha<-wt.alpha/sum(wt.alpha)
			alpha0<-sum(alpha.try*wt.alpha)
			#alpha0<-alpha.try[wt.alpha==max(wt.alpha)][1]
			#alpha0<-1

		}
		


		lambda.shrink<-lambda.all<-rep(grand.lambda,p)*wts.lam

		
		





	if(type=="probit"){
		
		loss.int<-function(int){
			(mean(pnorm(fits+int))-mean(y0))^2
		}
		
		
		probit.int<- optimize(f=loss.int,interval=c(-5,5))$min
		
		if(!EM){
		y[y0==1]<-rtnorm(n=sum(y0==1),mean=fits[y0==1]+probit.int,sd=1,lower=0,upper=Inf)
		y[y0==0]<-rtnorm(n=sum(y0==0),mean=fits[y0==0]+probit.int,sd=1,lower=-Inf,upper=0)
		}else{

		fits<-fits+probit.int
		alpha.tn<-0-fits[y0==1]
		beta.tn<-Inf
		num<-dnorm(beta.tn)-dnorm(alpha.tn)
		denom<-pnorm(beta.tn)-pnorm(alpha.tn)
		y[y0==1]<-fits[y0==1]-(num/denom)
		
		alpha.tn<--Inf
		beta.tn<-fits[y0==0]-0
		num<-dnorm(beta.tn)-dnorm(alpha.tn)
		denom<-pnorm(-beta.tn)-pnorm(alpha.tn)
		y[y0==0]<-fits[y0==0]-(num/denom)

		#probit.int<-mean(y)
		#y<-y-probit.int
		#y[fits>10]<-10
		#y[fits<-10]<- -10
		#y[!is.finite(y)]<-sign(fits[!is.finite(y)])*20

		


		y<-y-mean(y)
	#print(beta.curr)
			}
			sigma.sq<-1
	}
	
	if(type=="tobit"){
		if(!EM){
		fits.trunc<-X.big%*%beta.curr
		if(lower.trunc==TRUE) y[y0==trunc]<-rtnorm(sum(y0==trunc),lower=-Inf,upper=trunc,mean=fits.trunc[y0==trunc],sd=sigma.sq^.5)
		if(lower.trunc==FALSE) y[y0==trunc]<-rtnorm(sum(y0==trunc),lower=trunc,upper=Inf,mean=fits.trunc[y0==trunc],sd=sigma.sq^.5)
	y<-(y-mean(y))
		}
	if(EM){
		fits<-X.big%*%beta.curr
		int.adj<-mean(y[y0>trunc])-mean(y0[y0>trunc])
		y<-y-int.adj
		fits<-fits-int.adj
		#fits<-fits-mean(fits[y0>trunc])+mean(y[y0>trunc])
	
		if(lower.trunc==TRUE){
		alpha.tn<--Inf
		beta.tn<- -(y[y0==trunc][1]-fits[y0==trunc])/sigma.sq^.5#fits[y0==0]-0
		num<-dnorm(beta.tn)-dnorm(alpha.tn)
		denom<-pnorm(-beta.tn)-pnorm(alpha.tn)
		y[y0==trunc]<-fits[y0==trunc]-(num/denom)
		#print(y[y0==trunc])
		}
		y<-(y-mean(y))
		fits<-fits-mean(fits)
		
	}
	
	}
	}#End thin loop
	
	if(i.gibbs>burnin){
		beta.mean[,i.gibbs-burnin]<-beta.curr
		beta.mode[,i.gibbs-burnin]<-beta.curr.mode
	sigma.sq.run[i.gibbs-burnin]<-sigma.sq
	lambda.run[,i.gibbs-burnin]<-lambda.all
	beta.ci[,i.gibbs-burnin]<-beta.curr.ci

	}
	
	if(i.gibbs%in%seq(1,burnin+gibbs,length=8)&save.temp){
		save(file="temp_sparsereg",beta.mode)
	}
	
	
	if(!EM){
	if(i.gibbs==1) cat("    0% of Burnin Completed:", i.gibbs,"/",burnin,"\n")	
	if(i.gibbs==floor(burnin/4)) cat("   25% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(burnin/2)) cat("   50% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(3*burnin/4)) cat("   75% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(burnin)) cat("  100% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs == burnin+1) cat("Step 3 of 4: Burnin Completed; Saving Samples \n")
	if(i.gibbs==floor(burnin+1)) cat("    0% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	if(i.gibbs==floor(burnin+(gibbs)/4)) cat("   25% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	if(i.gibbs==floor(burnin+(gibbs)/2)) cat("   50% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	if(i.gibbs==floor(burnin+3*(gibbs)/4)) cat("   75% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	if(i.gibbs==floor(burnin+gibbs)) cat("  100% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	}
	##Gather things
	}#End gibbs loop

if(EM){
	if(i.gibbs%%10==0) cat("  EM iteration ", i.gibbs, "\n")
}

if(EM) cat("Step 3 of 3: Gathering Output")
if(!EM)	cat("Step 4 of 4: Gathering Output\n")
	
	if(type=="probit") scale.y<-1

		for(i.scale in 1:nrow(beta.mean)) {
			beta.mean[i.scale,]<-beta.mean[i.scale,]/scale.back[i.scale]*scale.y
			beta.mode[i.scale,]<-beta.mode[i.scale,]/scale.back[i.scale]*scale.y
			beta.ci[i.scale,]<-beta.ci[i.scale,]/scale.back[i.scale]*scale.y
			}
			

		
	
	matrix.names<-matrix(NA,nrow=3,ncol=ncol(X.big))
	rownames(matrix.names)<-c("X","treat1","treat2")

	for(i.scale in 1:ncol(X.big)) X.big[,i.scale]<-X.big[,i.scale]*scale.back[i.scale]


	row.names(beta.mode)<-gsub("_",": ",row.names(beta.mode))
	row.names(beta.mean)<-gsub("_",": ",row.names(beta.mean))

	if(EM){
		max.rows<-max(which(!is.na(beta.mean[1,])))
		beta.mean<-beta.mean[,max.rows]
		beta.mode<-beta.mode[,max.rows]
		beta.ci<-beta.mode
		beta.mean<-rbind(beta.mean,beta.mean,beta.mean)
		beta.mode<-rbind(beta.mode,beta.mode,beta.mode)
		beta.ci<-beta.mode
	}  else{
		beta.mode=as.mcmc(t(beta.mode))
		beta.mean=as.mcmc(t(beta.mean))
		beta.ci=as.mcmc(t(beta.ci))
		}
	

	output<-list("beta.mode"=beta.mode,"beta.mean"=beta.mean,"beta.ci"=beta.ci,
	"sigma.sq"=sigma.sq.run*scale.y^2,"X"=X.big,"y"=y,"lambda"=lambda.run,
	"varmat"=matrix.convert,"baseline"=baseline.vec,"modeltype"="onestage")
	class(output)<-"sparsereg"
	return(output)
}#End function loop
	
#sparsereg<-function(y, X,...) 	subgroup(y, X, treat=NULL, ...)	
#csts<-function(y, X, id=NULL, id2=NULL, id3=NULL, ...) 	subgroup(y, X,id=id,id2=id2,id3=id3,...)


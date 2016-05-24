glmm <-
function(fixed,random, varcomps.names,data, family.glmm, m,varcomps.equal, doPQL=TRUE, debug=FALSE,p1=1/3,p2=1/3,p3=1/3,rmax=1000,iterlim=1000,par.init=NULL,zeta=5){
	if(missing(varcomps.names)) stop("Names for the variance components must be supplied through varcomps.names")
	if(is.vector(varcomps.names)!=1) stop("varcomps.names must be a vector")
	if(missing(varcomps.equal)){
		varcomps.equal<- c(1:length(varcomps.names))}
	call<-match.call()

	#this much will figure out how to interpret the formula
	#first the fixed effects part
	stopifnot(inherits(fixed, "formula"))
    if (missing(data)) {
        barf <- lm(fixed, method = "model.frame")
    } else {
        stopifnot(inherits(data, "data.frame"))
        barf <- lm(fixed, data = data, method = "model.frame")
    }
    x <- model.matrix(fixed, data = barf)
    y <- model.response(barf)
    
    #then the part for the random effects. 
    #first, if it's not a list, make it a list 
    randcall<-random
	if (! is.list(random))
        random <- list(random)
    #put this stuff in a loop and loop along the list
    #for i in 1:length(formula2)
    for (irandom in seq(along = random)) {
    		r<-random[[irandom]]
    		stopifnot(inherits(r, "formula"))
    		if (missing(data)) {
       		 barf2 <- lm(r, method = "model.frame")
   		} else {
       		 stopifnot(inherits(data, "data.frame"))
        	barf2 <- lm(r, data = data, method = "model.frame")
    		}
    	random[[irandom]] <- model.matrix(r, data = barf2)
	#thisgroup<-varcomps.equal[irandom]
	#names(random)[irandom]<-varcomps.names[thisgroup]

	if(length(y)!=nrow(random[[irandom]])) stop("Fixed and random effect model matrices should have same number of rows")
	}
	#so now random is a list containing a model matrix for each formula, and some matrices share variance components

	if(is.numeric(varcomps.equal)==F) stop("varcomps.equal must be a vector containing numbers to indicate which variance components are equal.")
	if(length(varcomps.equal)!=length(random)){
		stop("The length of varcomps.equal must be equal to the length of the random-effects call.")} 
	if(length(unique(varcomps.equal))!=length(varcomps.names)){
		stop("You must name each unique variance component. Check varcomps.names and varcomps.equal.")} 
	if(min(varcomps.equal)!=1)stop("The vector varcomps.equal must contain numbers starting at 1 to denote which variance components are equal.")	
	levs<-ordered(unique(varcomps.equal))
	family.glmm<-getFamily(family.glmm)
	family.glmm$checkData(y)


	#check p1 p2 p3
	if(!is.numeric(p1))stop("p1 must be a number between 0 and 1")
	if(p1>1) stop("p1 must be a number between 0 and 1")
	if(p1<0) stop("p1 must be a number between 0 and 1")
	if(p1==0) stop("p1 must be nonzero")
	if(!is.numeric(p2))stop("p2 must be a number between 0 and 1")
	if(p2>1) stop("p2 must be a number between 0 and 1")
	if(p2<0) stop("p2 must be a number between 0 and 1")
	if(!is.numeric(p3))stop("p3 must be a number between 0 and 1")
	if(p3>1) stop("p3 must be a number between 0 and 1")
	if(p3<0) stop("p3 must be a number between 0 and 1")
	if(p1+p2+p3!=1) stop("p1+p2+p3 must equal 1")
	
	#this loop is a 2-4-1. We want to check that they're filling in varcomps.equal correctly. 
	#We also want to group all the design matrices that share a variance components.
	#Now z is a list with the number of design mats = number of distinct variance components
	z<-list()
	for(i in 1:length(levs)){
		if(levs[i]!=i) stop("The numbers in the vector varcomps.equal must be consecutive. You must start at 1 and then each entry must be the next consecutive number or a repeat of a previous number.")
		these<-varcomps.equal==i
		thesemats<-random[these]
		z[[i]]<-do.call(cbind,thesemats)
	}
	names(z)<-varcomps.names

	mod.mcml<-list(x = x, z=z,y = y)

	
	#so now the 3 items are x (matrix), z (list), y (vector)
	#end figuring out how to interpret the formula

	#make sure par.init has the right number of parameters and that the variance 
	#components are positive to start. also skip pql bc using par.init instead
	if(!is.null(par.init)){
		nbeta<-ncol(x)
		nbetaplusT<-nbeta+length(z)
		if(length(par.init)!=nbetaplusT) stop("par.init is not the correct length. It should contain initial values for the fixed effects and variance components.")
		vcs<-par.init[-(1:nbeta)]
		if(any(vcs<=10^-9)) stop("Initial values for the variance components in par.init must be positive and sufficiently large (greater than 10^-9).")
		doPQL<-FALSE
		#if par.init is given, we want to use those not PQL
	}
	
	#cache will hold some pql estimates and the importance sampling weights that wouldn't otherwise be returned
	cache <- new.env(parent = emptyenv())

	#if the user wants to do pql, do it and use that as the trust start point
	if(doPQL==TRUE){
	      #do PQL
	      pql.out<-pql(mod.mcml,family.glmm,cache)
	      s.pql<-cache$s.twid	
	      sigma.pql<-pql.out$sigma
	      nu.pql<-sigma.pql^2
	      beta.pql<-cache$beta.twid 
		  par.init<-c(beta.pql,nu.pql) 
	}
	
	#if the user does not want to do pql, then the best guess of the rand effs is 0
	#if the user did not provide par.init, then an arbitrary guess is 0 for beta
	# and 1 for nu
	if(doPQL==FALSE){
	    nrand<-lapply(mod.mcml$z,ncol)
	    nrandom<-unlist(nrand)
	    totnrandom<-sum(nrandom)
	    s.pql<-rep(0,totnrandom)
		if(!is.null(par.init)){ #then par.init is already specified by user
			beta.pql<-par.init[1:nbeta]
			nu.pql<-par.init[-(1:nbeta)]
			sigma.pql<-sqrt(nu.pql)
		}
		if(is.null(par.init)){
	        sigma.pql<-nu.pql<-rep(1,length(mod.mcml$z))
	        beta.pql<-rep(0,ncol(mod.mcml$x))
		    par.init<-c(beta.pql,nu.pql) 
		}
	}

	#calculate A*, D* and u*
	nrand<-lapply(mod.mcml$z,ncol)
	nrandom<-unlist(nrand)
	q<-sum(nrandom)
	if(q!=length(s.pql)) stop("Can't happen. Number of random effects returned by PQL must match number of random effects specified by model.")
	eek<-getEk(mod.mcml$z)
	#if any of the variance components are too close to 0, make them bigger:
	if(any(sigma.pql<10^-3)){
		theseguys<-which(sigma.pql<10^-3)
		sigma.pql[theseguys]<-10^-3
	}
	Aks<-Map("*",eek,sigma.pql)
	A.star<-addVecs(Aks) #at this point still a vector
	D.star<-A.star*A.star #still a vector
	u.star<-A.star*s.pql 
	Dstarinvdiag<-1/D.star
	Dstarnotsparse<-diag(D.star)
	D.star.inv<-Diagonal(length(u.star),Dstarinvdiag)
	D.star<-Diagonal(length(u.star),D.star)

	#now D.star.inv and D.star are both diagonal matrices
	#Diagonal from Matrix package is used bc these are sparse matrices
	#If q (# rand effs) is large, then need to be careful with these

	#determine m1, m2, m3 based on probs p1, p2, p3
	foo<-runif(m)
	m1<-sum(foo<p1)
	m2<-sum(foo<p1+p2)-m1	
	m3<-m-m1-m2

#	#generate m1 from N(0,I) and will be scaled to N(0,D) later
#	zeros<-rep(0,length(u.star))
#	ones<-rep(1,length(u.star))
#	ident<-Diagonal(length(u.star),ones)
#	genData<-genRand(zeros,ident,m1)
#	

	#generate m1 from t(0,D*)
	if(m1>0) genData<-rmvt(m1,sigma=Dstarnotsparse,df=zeta,type=c("shifted"))
	if(m1==0) genData<-NULL		

	#generate m2 from N(u*,D*)
	if(m2>0) genData2<-genRand(u.star,D.star,m2)
	if(m2==0) genData2<-NULL


	#generate m3 from N(u*,(Z'c''(Xbeta*+zu*)Z+D*^{-1})^-1)
	if(m3>0){
		Z=do.call(cbind,mod.mcml$z)
		eta.star<-as.vector(mod.mcml$x%*%beta.pql+Z%*%u.star)
		cdouble<-family.glmm$cpp(eta.star) #still a vector
		cdouble<-Diagonal(length(cdouble),cdouble)
		Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
		Sigmuh<-solve(Sigmuh.inv)
		genData3<-genRand(u.star,Sigmuh,m3)
	}
	if(m3==0) genData3<-NULL

#	#these are from distribution based on data
#	if(distrib=="tee")genData<-genRand(sigma.gen,s.pql,mod.mcml$z,m1,distrib="tee",gamm)
#	if(distrib=="normal")genData<-genRand(sigma.pql,s.pql,mod.mcml$z,m1,distrib="normal",gamm)
#	#these are from standard normal
#	ones<-rep(1,length(sigma.pql))
#	zeros<-rep(0,length(s.pql))
#	genData2<-genRand(ones,zeros,mod.mcml$z,m2,distrib="normal",gamm)

	umat<-rbind(genData,genData2,genData3)

	#use trust to max the objfun (monte carlo likelihood)
	trust.out<-trust(objfun,parinit=par.init,rinit=10, minimize=FALSE, rmax=rmax, iterlim=iterlim, blather=debug, nbeta=length(beta.pql), nu.pql=nu.pql, 
umat=umat, mod.mcml=mod.mcml, family.glmm=family.glmm, u.star=u.star,  cache=cache,  p1=p1,p2=p2, p3=p3,m1=m1, D.star=D.star,Sigmuh=Sigmuh,Sigmuh.inv=Sigmuh.inv,zeta=zeta)

	beta.trust<-trust.out$argument[1:length(beta.pql)]
	nu.trust<-trust.out$argument[-(1:length(beta.pql))]

#while(trust.out$converged==FALSE){

#}

	names(beta.trust)<-colnames(mod.mcml$x)
	names(nu.trust)<-varcomps.names

	if(debug==TRUE){
	debug<-list(beta.pql=beta.pql, nu.pql=nu.pql,  trust.argpath=trust.out$argpath, u.star=u.star, umat=umat,weights=cache$weights,wtsnumer=cache$numer,wtsdenom=cache$denom,m1=m1,m2=m2,m3=m3,trust.argtry=trust.out$argtry, trust.steptype=trust.out$steptype, trust.accept=trust.out$accept, trust.r=trust.out$r, trust.rho=trust.out$rho, trust.valpath=trust.out$valpath, trust.valtry=trust.out$valtry, trust.preddif=trust.out$preddif, trust.stepnorm=trust.out$stepnorm)
	}
	
	return(structure(list(beta=beta.trust,nu=nu.trust, likelihood.value=trust.out$value, likelihood.gradient=trust.out$gradient, likelihood.hessian=trust.out$hessian,
	trust.converged=trust.out$converged,  mod.mcml=mod.mcml,
	 fixedcall=fixed,randcall=randcall, x=x,y=y, z=random,
	family.glmm=family.glmm, call=call, varcomps.names=varcomps.names, 
	varcomps.equal=varcomps.equal, debug=debug), class="glmm"))
}

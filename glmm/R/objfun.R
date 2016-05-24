objfun <-
function(par,nbeta,nu.pql,umat,u.star,mod.mcml,family.glmm,cache,p1,p2,p3,m1,D.star, Sigmuh, Sigmuh.inv, zeta){

	beta<-par[1:nbeta]
	nu<-par[-(1:nbeta)]
	m<-nrow(umat)

	if (!missing(cache)) stopifnot(is.environment(cache))

	if(any(nu<=0)){
		out<-list(value=-Inf,gradient=rep(1,length(par)),hessian=as.matrix(c(rep(1,length(par)^2)),nrow=length(par)))
	return(out)
	}
	
	Z=do.call(cbind,mod.mcml$z)
	T<-length(mod.mcml$z)
	nrand<-lapply(mod.mcml$z,ncol)
	nrandom<-unlist(nrand)


	family.glmm<-getFamily(family.glmm)
	if(family.glmm$family.glmm=="bernoulli.glmm"){family_glmm=1}	
	if(family.glmm$family.glmm=="poisson.glmm"){family_glmm=2}	

	Dstarinvdiag<-1/diag(D.star)
	D.star.inv<-diag(Dstarinvdiag)

	logdet.D.star.inv<-	-sum(log(diag(D.star)))
	logdet.Sigmuh.inv<-sum(log(eigen(Sigmuh.inv,symmetric=TRUE)$values))
 	myq<-nrow(D.star.inv)

	tconst<-tconstant(zeta,myq,Dstarinvdiag)

	#for the particular value of nu we're interested in, need to prep for distRandGenC
	eek<-getEk(mod.mcml$z)
	preDinvfornu<-Map("*",eek,(1/nu))
	Dinvfornu<-addVecs(preDinvfornu)
	logdetDinvfornu<-sum(log(Dinvfornu))
	Dinvfornu<-diag(Dinvfornu)
	
	meow<-rep(1,T+1)
	meow[1]<-0
	throwaway<-T+1
	meow[2:throwaway]<-cumsum(nrandom)
	
	pea<-c(p1,p2,p3)
	n<-nrow(mod.mcml$x)

##need to scale first m1 vectors of generated random effects by multiplying by A

#	preAfornu<-Map("*",eek,sqrt(nu))
#	Afornu<-addVecs(preAfornu)

#	for(k in 1:m1){
#		u.swoop<-umat[k,]
#		umat[k,]<-u.swoop*Afornu
#		}

	stuff<-.C("objfunc", as.double(mod.mcml$y),as.double(t(umat)), as.integer(myq), as.integer(m), as.double(mod.mcml$x), as.integer(n), as.integer(nbeta), as.double(beta), as.double(Z), as.double(Dinvfornu), as.double(logdetDinvfornu),as.integer(family_glmm), as.double(D.star.inv), as.double(logdet.D.star.inv), as.double(u.star), as.double(Sigmuh.inv), as.double(logdet.Sigmuh.inv), pea=as.double(pea), nps=as.integer(length(pea)), T=as.integer(T), nrandom=as.integer(nrandom), meow=as.integer(meow),nu=as.double(nu), zeta=as.integer(zeta),tconst=as.double(tconst), v=double(m),value=double(1),gradient=double(length(par)),hessian=double((length(par))^2),PACKAGE="glmm")


	if (!missing(cache)) cache$weights<-stuff$v		

	list(value=stuff$value,gradient=stuff$gradient,hessian=matrix(stuff$hessian,ncol=length(par),byrow=FALSE))

}


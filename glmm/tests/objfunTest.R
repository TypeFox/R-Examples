#check objfun using finite differences
library(glmm)
data(BoothHobert)
set.seed(1234)
out<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE)
mod.mcml<-out$mod.mcml
debug<-out$debug
nu.pql<-debug$nu.pql
nu.pql
beta.pql<-debug$beta.pql
beta.pql
family.glmm<-out$family.glmm
umat<-debug$umat
u.pql<-debug$u.star
m1<-debug$m1

par<-c(6,1.5)
del<-rep(10^-8,2)

objfun<-glmm:::objfun
getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

#need to get D*
eek<-getEk(mod.mcml$z)
Aks<-Map("*",eek,nu.pql)
D.star<-addVecs(Aks) 
D.star<-diag(D.star)
D.star.inv<-solve(D.star)


#need to also recreate the variance matrix of last imp sampling distribution
Z=do.call(cbind,mod.mcml$z)
eta.star<-as.vector(mod.mcml$x%*%beta.pql+Z%*%u.pql)
cdouble<-bernoulli.glmm()$cpp(eta.star) #still a vector
cdouble<-diag(cdouble)
Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
Sigmuh<-solve(Sigmuh.inv)

p1=p2=p3=1/3
zeta=5

# define a few things that will be used for finite differences
lth<-objfun(par=par, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,p1=p1,p2=p2,p3=p3,m1=m1, Sigmuh=Sigmuh, D.star=D.star, Sigmuh.inv= Sigmuh.inv, zeta=zeta)
lthdel<-objfun(par=par+del, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,p1=p1,p2=p2,p3=p3,m1=m1, Sigmuh=Sigmuh, D.star=D.star, Sigmuh.inv=Sigmuh.inv, zeta=zeta)

all.equal(as.vector(lth$gradient%*%del),lthdel$value-lth$value)
all.equal(as.vector(lth$hessian%*%del),lthdel$gradient-lth$gradient)

#see exactly how big the difference is
as.vector(lth$gradient%*%del)-(lthdel$value-lth$value)
as.vector(lth$hessian%*%del)-(lthdel$gradient-lth$gradient)

#we know these differences are small when we compare it to the actual values
lthdel$value-lth$value
as.vector(lth$gradient%*%del)
as.vector(lth$hessian%*%del)
lthdel$gradient-lth$gradient

##########################################
##### to make sure that the objfun function is correct, compare it against the version without any C code. here is objfun without c:
objfunNOC <-
function(par,nbeta, nu.pql,umat, u.star=u.star, mod.mcml,family.glmm, cache,gamm,p1,p2,p3, D.star, Sigmuh, zeta){

	#print(par)
	beta<-par[1:nbeta]
	nu<-par[-(1:nbeta)]
	D<-nu*diag(10)
	D.inv<-(1/nu)*diag(10)
	m<-nrow(umat)

	if (!missing(cache)) stopifnot(is.environment(cache))
	if(missing(cache)) cache<-new.env(parent = emptyenv())

	if(sum(nu<=0)>0){
		out<-list(value=-Inf,gradient=rep(1,length(par)),hessian=as.matrix(c(rep(1,length(par)^2)),nrow=length(par)))
	return(out)
	}
	
	Z=do.call(cbind,mod.mcml$z)

	eta<-b<-rep(0,m)
	lfu<-lfyu<-list(rep(c(0,0,0),m))
	lfu.twid<-matrix(data=NA,nrow=m,ncol=4)

	D.star.inv<-solve(D.star)
	Sigmuh.inv<-solve(Sigmuh)
	Dstinvdiag<-diag(D.star.inv)
	tconst<-tconstant(zeta,nrow(D.star.inv),Dstinvdiag)
	
	#for each simulated random effect vector
	for(k in 1:m){
		Uk<-umat[k,]  #use the simulated vector as our random effect vec
		eta<-mod.mcml$x%*%beta+Z%*%Uk # calculate eta using it
		zeros<-rep(0,length(Uk))

		#log f_theta(u_k)
		lfu[[k]]<-distRand(nu,Uk,mod.mcml$z,zeros) 

		#log f_theta(y|u_k)
		lfyu[[k]]<-elR(mod.mcml$y,mod.mcml$x,eta,family.glmm) 

		#log f~_theta(u_k)
		lfu.twid[k,1]<-tdist2(tconst,Uk,D.star.inv,zeta=zeta,myq=nrow(D.star.inv))
		lfu.twid[k,2]<-distRandGeneral(Uk,u.star,D.star.inv)
		lfu.twid[k,3]<-distRandGeneral(Uk,u.star,Sigmuh.inv)
		
		tempmax<-max(lfu.twid[k,1:3])
		blah<-exp(lfu.twid[k,1:3]-tempmax)
		pea<-c(p1,p2,p3)
		qux<-pea%*%blah
		lfu.twid[k,4]<-tempmax+log(qux)
		
		b[k]<-as.numeric(lfu[[k]]$value)+as.numeric(lfyu[[k]]$value)-lfu.twid[k,4]
	}

	a<-max(b)
	thing<-exp(b-a)
	value<-a-log(m)+log(sum(thing))
	v<-thing/sum(thing)
	#bk are log weights
	cache$weights<-exp(b)
	
	Gpiece<-matrix(data=NA,nrow=nrow(umat),ncol=length(par))
	
	#lfuky<-NA
	for(k in 1:nrow(umat)){
		Gpiece[k,]<-c(lfyu[[k]]$gradient,lfu[[k]]$gradient)*v[k]	
		
				#lfuky[k]<-c(lfyu[[k]]$gradient,lfu[[k]]$gradient)
		#Gpiece[k,]<-lfuky[k]*v[k]	
	}
	G<-apply(Gpiece,2,sum)
	
	#Hessian has three pieces: panda, lobster, GGT
	panda.list<-list()
	for(k in 1:nrow(umat)){
		panda.list[[k]]<-c(lfyu[[k]]$gradient,lfu[[k]]$gradient)%*%t(c(lfyu[[k]]$gradient,lfu[[k]]$gradient))*v[[k]]

	}
	panda<-addMats(panda.list)
	
	lobster.list<-list()
	for(k in 1:nrow(umat)){

		mat1<-lfyu[[k]]$hessian
		mat2<-lfu[[k]]$hessian

		d1<-nrow(mat1)
		d2<-nrow(mat2)
		newmat<-matrix(data=0,nrow=d1+d2,ncol=d1+d2)

		newmat[1:d1,1:d1]<-mat1
		here<-d1+1
		there<-d1+d2
		newmat[here:there,here:there]<-mat2	
		lobster.list[[k]]<-newmat*v[k]
	}
	lobster<-addMats(lobster.list)
	
	hessian<-lobster+panda-G%*%t(G)
	list(value=value,gradient=G,hessian=hessian)
}

#here is el without C
elR <-
function(Y,X,eta,family.mcml){
	family.mcml<-getFamily(family.mcml)
	neta<-length(eta)

	if(family.mcml$family.glmm=="bernoulli.glmm"){
		foo<-.C("cum3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cumout=double(1))$cumout
		mu<-.C("cp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cpout=double(neta))$cpout
		cdub<-.C("cpp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cppout=double(neta))$cppout
	}
	if(family.mcml$family.glmm=="poisson.glmm"){
		foo<-.C("cum3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cumout=double(1))$cumout
		mu<-.C("cp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cpout=double(neta))$cpout
		cdub<-.C("cpp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cppout=double(neta))$cppout
	}

	value<-as.numeric(Y%*%eta-foo)
	gradient<-t(X)%*%(Y-mu)	
	cdubmat<-diag(cdub)
	hessian<-t(X)%*%(-cdubmat)%*%X
	
	list(value=value,gradient=gradient,hessian=hessian)
}
#here are some other functions we'll need to compare objfun and objfunNOC
getFamily<-glmm:::getFamily
addMats<-glmm:::addMats
tdist2<-function(tconst,u, Dstarinv,zeta,myq){
	inside<-1+t(u)%*%Dstarinv%*%u/zeta
	logft<-tconst - ((zeta+myq)/2)*log(inside)
	as.vector(logft)
}

tconstant<-glmm:::tconstant
distRandGeneral<-function(uvec,mu,Sigma.inv){
	logDetSigmaInv<-sum(log(eigen(Sigma.inv,symmetric=TRUE)$values))
	umu<-uvec-mu
	piece2<-t(umu)%*%Sigma.inv%*%umu
	out<-as.vector(.5*(logDetSigmaInv-piece2))
	const<-length(uvec)*.5*log(2*pi)
	out<-out-const
	out
}


distRand <-
function(nu,U,z.list,mu){
	# T=number variance components
	T<-length(z.list)
	
	#nrandom is q_t
	nrand<-lapply(z.list,ncol)
	nrandom<-unlist(nrand)
	totnrandom<-sum(nrandom)
	
	mu.list<-U.list<-NULL
	if(T==1) {
		U.list[[1]]<-U
		mu.list[[1]]<-mu
		}

	if(T>1){
		U.list[[1]]<-U[1:nrandom[1]] 
		mu.list[[1]]<-mu[1:nrandom[1]]
		for(t in 2:T){
			thing1<-sum(nrandom[1:t-1])+1
			thing2<-sum(nrandom[1:t])
			U.list[[t]]<-U[thing1:thing2]
			mu.list[[t]]<-mu[thing1:thing2]
		}
	}
	
	val<-gradient<-Hessian<-rep(0,T)
	
	#for each variance component
	for(t in 1:T){
		you<-as.vector(U.list[[t]])
		mew<-as.vector(mu.list[[t]])
		Umu<-(you-mew)%*%(you-mew)
		val[t]<- -length(U)*log(2*pi)/2+as.numeric(-.5*nrandom[t]*log(nu[t])-Umu/(2*nu[t]))
		
		gradient[t]<- -nrandom[t]/(2*nu[t])+Umu/(2*(nu[t])^2)
		
		Hessian[t]<- nrandom[t]/(2*(nu[t])^2)- Umu/((nu[t])^3)
		
	}
		
	value<-sum(val)
	if(T>1) hessian<-diag(Hessian)
	if(T==1) hessian<-matrix(Hessian,nrow=1,ncol=1)
	
	list(value=value,gradient=gradient,hessian=hessian)		
}

#finally, compare objfun and objfunNOC for B+H example

that<-objfunNOC(par=par, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,p1=p1,p2=p2,p3=p3, Sigmuh=Sigmuh,D.star=D.star, zeta=zeta)
all.equal(that,lth)

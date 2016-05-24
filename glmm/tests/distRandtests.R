library(glmm)

nrandom<-c(4,2,2)
T<-length(nrandom)

#create meow
meow<-rep(0,T+1)
	meow[1]=0
	meow[2]=nrandom[1]
	if(T>2){	
		for(t in 2:T+1){
		meow[t]=meow[t-1]+nrandom[t-1]
		}
	}

set.seed(1234)
U<-rnorm(8)
nu<-c(.5,.3,.9)
mu<-rep(0,8)

#third version of distRand used in objfun
#calculates only gradient and hessian 
#value is done in distRandGenC
#distRandGenC is tested in testpiecesBH.R
distRand3 <-
function(nu,U,mu,T,nrandom,meow){
	# T=number variance components
	#nrandom is q_t
	#meow will help us figure out which t working with

	Umu<-gradient<-hessvec<-rep(0,T)
	
	#for each variance component
	for(t in 1:T){
		#need to calculate (Ut-mut)'(Ut-mut)
		for(i in (meow[t]+1):meow[t+1]){
			Umu[t]=Umu[t]+(U[i]-mu[t])^2
		}
		gradient[t]=Umu[t]/(2*(nu[t])^2)-nrandom[t]/(2*nu[t])
		hessvec[t]=	-Umu[t]/((nu[t])^3)+nrandom[t]/(2*(nu[t])^2)
	}
		

	if(T>1) hessian<-diag(hessvec)
	if(T==1) hessian<-matrix(hessvec,nrow=1,ncol=1)
	
	list(gradient=gradient,hessian=hessian)		
}

#second version of distRand
distRand2 <-
function(nu,U,mu,T,nrandom){
	# T=number variance components
	#nrandom is q_t
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
		val[t]<- as.numeric(-.5*nrandom[t]*log(nu[t])-Umu/(2*nu[t]))
		
		gradient[t]<- -nrandom[t]/(2*nu[t])+Umu/(2*(nu[t])^2)
		
		Hessian[t]<- nrandom[t]/(2*(nu[t])^2)- Umu/((nu[t])^3)
		
	}
		
	value<-sum(val)
	if(T>1) hessian<-diag(Hessian)
	if(T==1) hessian<-matrix(Hessian,nrow=1,ncol=1)
	
	list(value=value,gradient=gradient,hessian=hessian)		
}



dr3<-distRand3(nu,U,mu,T,nrandom,meow)
dr2<-distRand2(nu,U,mu,T,nrandom)

all.equal(dr3$gradient,dr2$gradient)
all.equal(dr3$hessian,dr2$hessian)

drc<-.C("distRand3C", as.double(nu), as.double(mu), as.integer(T), as.integer(nrandom), as.integer(meow), as.double(U), double(T), double(T^2))

all.equal(drc[[7]],dr3$gradient)
all.equal(matrix(drc[[8]],nrow=3,byrow=F),dr3$hessian)

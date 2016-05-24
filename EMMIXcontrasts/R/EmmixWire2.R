### start



#E-step

tau.estep.wire<-function(dat,pro,mu,sigma,n,m,g)
{		
obj <- .Fortran("estepmvn",PACKAGE="EMMIXcontrasts",
	as.double(dat),as.integer(n),as.integer(m),as.integer(g),
	as.double(pro),as.double(mu),as.double(sigma),
	mtauk=double(n*g),double(g),loglik=double(1),
	error = integer(1))[8:11] 
	if(obj$error) stop("error")
	tau <- array(obj$mtauk,c(n,g))
list(tau=tau,loglik=obj$loglik,pro=colMeans(tau))		
}

wire.init.fit<-function(dat,X,qe,n,m,g,nkmeans,nrandom=0)
{
wire.init.km<-function(dat,X,qe,n,m,g)
{
cluster<-rep(1,n)		
if(g>1)
cluster<- kmeans(dat,g,nstart=5)$cluster
wire.init.reg(dat,X,qe,n,m,g,cluster)
}	
wire.init.rd<-function(dat,X,qe,n,m,g)
{
cluster<-rep(1,n)		
if(g>1)
cluster<- sample(1:g,n,replace=TRUE)
wire.init.reg(dat,X,qe,n,m,g,cluster)
}
	
found<-NULL
found$loglik<--Inf
	
if(nkmeans>0) {
for(j in 1:nkmeans)
{	
initobj<-try(wire.init.km(dat,X,qe,n,m,g))		
if(length(initobj)>=4)
if(initobj$loglik>found$loglik)
found<-initobj			
} #end of loop
}
if(nrandom>0) {
for(j in 1:nrandom)
{
	initobj<-try(wire.init.rd(dat,X,qe,n,m,g))			
	if(length(initobj)>=4)
	if(initobj$loglik>found$loglik)
	found<-initobj
} #end of loop
}
found
}

wire.init.reg<-function(dat,X,qe,n,m,g,cluster)
{
	# set x for regression
	xx<-as.matrix(X)
	beta<-matrix(0,nrow=ncol(xx),ncol=g)
	sigma<-pro<-rep(0,g)
	mu     <- array(0,c(m,g))
	msigma <- array(0,c(m,m,g))
	lk <- rep(0,g)
	for( ij in 1:g)
	{
		ni<-sum(cluster==ij)
		if(ni==0){
			warning("empty cluster found!")
			next
		}
		nn     <- ni
		#pile up y
		y      <- c(t(dat[cluster==ij,]))	
		if(nn > 1)
		{
			for(i in 2:nn)
				xx <- rbind(xx,X)
		}
		obj<-lm(y~0+xx)
		beta[,ij]<-coef(obj)
		sigma[ij]<-sum((resid(obj))^2)/(nn*m)	
		pro[ij]<-ni/n
		xx<-as.matrix(X) # reset x for next iteration
		
	}
	
	# loglikelihood
	
	for(h in 1:g) {
		mu[,h] <-c(X%*%beta[,h])
		msigma[,,h] <- diag(sigma[h],m)
		ni <- sum(cluster==h)
		if(ni==0){
			warning("empty cluster found!")
			next
		}
	}
        ooo <- tau.estep.wire(dat,pro,mu,msigma,n,m,g)
	loglik <- ooo$loglik
	sigma.e<-matrix(0,ncol=g,nrow=qe)
	for(i in 1:qe)
		sigma.e[i,]<-sigma
	list(beta=beta,sigma.e=sigma.e,pro=pro,loglik=loglik,lk=lk)
}

emmixwire<-function(dat,g=1,ncov=3,nvcov=1,n1=0,n2=0,n3=0,
X=NULL,W=NULL,U=NULL,V=NULL,
cluster=NULL,init=NULL,debug=0,itmax=1000,epsilon=1e-5,nkmeans=5,nrandom=0)
{

mat.ABC.wire<-function(U,VV,W,sigma.e,DU,sigma.c,g,m)
{
		A<-B<-C<-BC<-array(0,dim=c(m,m,g))
		for(h in 1:g)
		{
			if(ncol(W)>1)
				A[,,h]<-diag(as.vector(W%*%sigma.e[,h]))
			else
			{
				A[,,h]<-diag(c(W[,1]*sigma.e[1,h]))
				
			}
			B[,,h]<-A[,,h]+U%*%DU[,,h]%*%t(U)
			C[,,h]<-sigma.c[h]*VV
			BC[,,h]<-B[,,h]+C[,,h]
		}
		list(A=A,B=B,C=C,BC=BC)
}
matM.wire<-function(B,C,tau,g,m)
{
		M<-array(0,dim=c(m,m,g))
		for(h in 1:g)
			M[,,h]<-solve(B[,,h]+C[,,h]*sum(tau[,h]))
		M
}
#--------------------------------------------
	
eb.estep<-function(tau,y,mu,DU,U,B,C,M,g,n,m,qb)
{
		eb1<-array(0,dim=c(n,qb,g))
		eb2<-array(0,c(g,qb,qb))
		invB<-array(0,c(m,m,g))
		
		for(h in 1:g)
		{
			invB[,,h]<-solve(B[,,h])
			# E(bi|y)
			
			eb1[,,h]<-t( DU[,,h]%*%t(U)%*%invB[,,h]%*%(  (t(y)-mu[,h])
																										-c(C[,,h]%*%M[,,h]%*%colSums(t(t(y)-mu[,h])*tau[,h])) ))
			#E(bi%*%bi^t|y)
			eb2[h,,]<-( sum(tau[,h])*DU[,,h]-(sum(tau[,h])-1)*DU[,,h]%*%t(U)%*%invB[,,h]%*%U%*%DU[,,h]	
									-DU[,,h]%*%t(U)%*%M[,,h]%*%U%*%DU[,,h])
	DU[,,h]<-(eb2[h,,]+ t(eb1[,,h]*tau[,h])%*%eb1[,,h])/sum(tau[,h])
}
		
list(DU=DU,eb1=eb1,invB=invB)
}
	
ec.estep<-function(tau,y,mu,sigma.c,V,M,g,n,qc)
{
		ec1<-array(0,c(qc,g))
		ec2<-kc<-rep(0,g)
		#
		for(h in 1:g)
		{
			ec1[,h]<-sigma.c[h]*t(V)%*%M[,,h]%*%colSums(t(t(y)-mu[,h])*tau[,h])
			kc[h] <-sigma.c[h]*qc-sum(diag(t(V)%*%M[,,h]%*%V))*(sigma.c[h]^2*sum(tau[,h]))
			ec2[h]<-c(t(ec1[,h])%*%ec1[,h]+kc[h])/qc
		}	
		#
		list(ec2=ec2,ec1=ec1)
	}
	
ee.estep<-function(y,mu,tau,U,V,W,A,invB,M,g,n,m,qe,dw,eb1,ec1)
{
ae<-ee<-array(0,dim=c(n,m,g))
ke<-rep(0,g)
thet<-matrix(0,ncol=g,nrow=qe)
mi<-diag(t(W)%*%W)
for(h in 1:g)
{
ae[,,h]<-t(mu[,h]+U%*%t(eb1[,,h])+c(V%*%ec1[,h]))			
ee[,,h]<- (y-ae[,,h])
for(id in 1:ncol(W))
{
AL<-A[,,h]*W[,id]
ke[h]<-(sum(tau[,h])*sum(diag(AL))-sum(diag(  t(AL)%*%M[,,h]%*%AL))
	-(sum(tau[,h])-1)*sum(diag( t(AL)%*%invB[,,h]%*%AL)))
		
	thet[id,h]<-( sum(c(  t((t(ee[,,h])*W[,id])^2)*tau[,h]  ))+ke[h])/(mi[id]*sum(tau[,h]))
	}# end of loop
			
	} #end of loop
	list(sigma.e=thet,ae=ae)
	}

	
tau2cluster<-function(tau)
{
apply(tau,FUN=which.max,MARGIN=1)
}

getcov <-function(msigma,sumtau,n,m,g,ncov)
{
sigma<-array(0,c(m,m))
if( (ncov==1)|(ncov==2))
{
for(h in 1:g)
sigma<-sigma+sumtau[h]*msigma[,,h]
sigma<-as.matrix(sigma/n)

if(ncov==2)
sigma<-diag(c(diag(sigma)),m)
for(h in 1:g)
msigma[,,h]=sigma
}

if(m>1)
{
if(ncov==4)
for(h in 1:g)
msigma[,,h]<-diag(c(diag(msigma[,,h])),m)

if(ncov==5)
for(h in 1:g)
msigma[,,h]<-diag(sum(diag(msigma[,,h]))/m,m)
}

msigma
}


# start EMMIX-WIRE analysis from initial values
fit.emmix.wire<-function(dat,X,W,U,V,pro,beta,sigma.e,sigma.b,sigma.c,
		 n,m,g,nb,qb,qc,qe,
		 debug,ncov,nvcov,itmax,epsilon,log=TRUE)
{
		# controls
		flag<-error<-0
		lk<-rep(0,itmax)
        
		oldpro<-pro
		
                nbeta<-beta

		
		VV<-V%*%t(V)
		dw<-diag(t(W)%*%W)
		xxx<-solve(t(X)%*%X)%*%t(X)

		mu<-matrix(0,ncol=g,nrow=m)
		
		for(i in 1:itmax)
		{
			mobj<-mat.ABC.wire(U,VV,W,sigma.e,sigma.b,sigma.c,g,m)
			A<-mobj$A
			B<-mobj$B
			C<-mobj$C
			BC<-mobj$BC
			for(h in 1:g) {
				mu[,h]     <- as.vector(X%*%beta[,h])
			}
			# E-step
			eobj<-tau.estep.wire(dat,oldpro,mu,BC,n,m,g)

			pro   <-eobj$pro
			tau   <-eobj$tau
			lk[i] <-eobj$loglik
			sumtau<- colSums(tau)
			M<-matM.wire(B,C,tau,g,m)
			obj1 <- eb.estep(tau,dat,mu,sigma.b,U,B,C,M,g,n,m,qb)
			obj2 <- ec.estep(tau,dat,mu,sigma.c,V,M,g,n,qc)
			obj3 <- ee.estep(dat,mu,tau,U,V,W,A,obj1$invB,M,g,n,m,qe,dw,obj1$eb1,obj2$ec1)
			
			
			# M-step
			
			if(ncov>0)
				sigma.b<-obj1$DU
			else
				for(h in 1:g)
					sigma.b[,,h]<-diag(0,qb)
			
			if( (ncov>0) & (ncov!=3) & (ncov!="AR") )
				sigma.b <- getcov(sigma.b,sumtau,n,qb,g,ncov)
			
			
			if(nvcov>0)
				sigma.c<-obj2$ec2
			else
				sigma.c<-rep(0,g)
			
			sigma.e<-obj3$sigma.e
			
			#--------------------------------
			
			for(h in 1:g)
				nbeta[,h]<-(beta[,h]+xxx%*%M[,,h]%*%A[,,h]%*%colSums(t(t(dat)-mu[,h])*tau[,h])/sumtau[h])
			
			#--------------------------------
			#
			loglik <- lk[i]

			if(debug)
			{
				cat('\n',i,'th,loglik=',loglik)
				
			}

beta   <- nbeta;
oldpro <- pro

if(i <= 10) next			
if(log){
if(abs(lk[i]-lk[i-10])<epsilon*abs(lk[i-10])) {flag<-1;break}
} else {
if(max(abs(c(nbeta-beta,pro-oldpro))) < epsilon) {flag<-1;break}
}

} # end of loop 
		
		if(flag==0) error <- 1
		
		#get the the final partition
		
		for(h in 1:g) {
			
			mu[,h]<-as.vector(X%*%beta[,h]+V%*%obj2$ec1[,h])
		}
		
		eobj2<-tau.estep.wire(dat,pro,mu,B,n,m,g)
		cluster<-tau2cluster(eobj2$tau)
	
		# 
# BIC & AIC
nu<-switch(paste(ncov),
	 '0' = 0,              # without u
	 '1' = qb*(1+qb)/2,    #common covariance
	 '2' = qb,             #common diagonal covariance
	 '3' = g*(qb*(1+qb)/2),#general covariance
	 '4' = g*qb,           #general diagonal covariance
	 '5' = g  )            #sigma(h)*I_m
		
nv<-ifelse(nvcov>0,g,0)
np<-(g-1)+nb*g + nu + nv+ g*qe
#
BIC<--2*loglik+np*log(n)
AIC<--2*loglik+np*2
		
if(debug)
cat('\n',g,"BIC=",BIC,"AIC=",AIC,
"\nloglik=",loglik,"np=",np)
		
		
#return values
		
ret <- list(error=error,loglik=loglik,np=np,
       BIC=BIC,AIC=AIC,cluster=cluster,
       pro=pro,beta=beta,sigma.e = sigma.e)
ret$lk =lk[lk!=0]
ret$tau=eobj2$tau
if(ncov==1||ncov==2||ncov==3||ncov==4||ncov==5)
{
ret$sigma.b <- sigma.b
ret$eb      <- obj1$eb1
}

if(nvcov>0)
{
ret$sigma.c <- sigma.c
ret$ec      <- obj2$ec1
}

#######################
ret
}
	
	
	
###############	

dat<-as.matrix(dat)
n<-nrow(dat)
m<-ncol(dat)

#
if(n1 >0 && n2>0 )
{
if(n3==0) {
X<-U<-cbind(rep(c(1,0),c(n1,n2)),rep(c(0,1),c(n1,n2)))
} else {
X<-U<-cbind(rep(c(1,0,0),c(n1,n2,n3)),rep(c(0,1,0),c(n1,n2,n3)),rep(c(0,0,1),c(n1,n2,n3)))
}
W <-rep(1,m)
V<-diag(m)
} else {

# check the matrix U
if(ncov==0)
{
U     <- diag(m)
} else {
if(is.null(U))
U <- cbind(rep(1,m))
}

# check the matrix V
if(is.null(V) || nvcov==0)
V  <- diag(m)

# check the matrix W
if(is.null(W))
	W <- cbind(rep(1,m))
}

#
# check the matrix X
if(is.null(X))
stop("X must be specified")

X<-as.matrix(X);U<-as.matrix(U)
V<-as.matrix(V);W<-as.matrix(W)
qb<-ncol(U);qc<-ncol(V);qe<-ncol(W);nb<-ncol(X)
##################################
#  some variables
#
tuv <- 0.2
# initialize the sigma_b and sigma_c
sigma.b<-array(0,c(qb,qb,g))
for(h in 1:g)
{
if(qb>1)
diag(sigma.b[,,h])<-rep(tuv,qb)
else
sigma.b[1,1,h] <- tuv
}
sigma.c<-rep(tuv,g)

#part 1: initial values
if(!is.null(init)){
	found=init
} else {
if(is.null(cluster))
{
found <- wire.init.fit(dat,X,qe,n,m,g,nkmeans,nrandom)
} else
found <- wire.init.reg(dat,X,qe,n,m,g,cluster)
}

if(length(found)<4)
	stop("not found inital values")

#part 2: call the main estimate procedure

ret<-fit.emmix.wire(dat,X,W,U,V,
found$pro,found$beta,found$sigma.e,sigma.b,sigma.c,
n,m,g,nb,qb,qc,qe,
debug,ncov,nvcov,itmax,epsilon)


if(qb==m && (ncov==4 || ncov==2))
{
tmp <- array(0,c(m,g))
for(h in 1:g)
tmp[,h] <- diag(ret$sigma.b[,,h])
ret$sigma.b <- tmp
}

if(ncov==5)
{
tmp <- rep(0,g)
for(h in 1:g)
tmp[h] <- diag(ret$sigma.b[,,h])[1]
ret$sigma.b <- tmp
}
ret$g=g
ret$m=m
ret$nb=nb
ret$X=X
ret$W=W
ret$U=U
ret$V=V

ret
}
######################


#differentially expressed genes (DEG)---   gene ranking

eq8.wire <-function(m,g,nb,X,W,U,V,sigma.e,sigma.b,sigma.c,nh,contrast)
{
omega <- rep(0,g)
if(is.null(contrast)) {
	if(nb==2) {
	  K1 = c(1,-1,rep(0,m))
	  K2 = c(1,-1)
	} else {
	  if(nb==3)
	  { 
		K1 = c(1,0,-1,rep(0,m))
		K2 = c(1,0,-1)
	  }
       }
} else {
        if(nb==2) {
          if(length(contrast)!=2)
           stop("contrast should be a vector with length of 2")
	  K1 = c(contrast,rep(0,m))
	  K2 = c(contrast)
        } else {
	  if(nb==3)
	  {
            if(length(contrast)!=3)
               stop("contrast should be a vector with length of 3")
            K1 = c(contrast,rep(0,m))
	    K2 = c(contrast)
	  }
       }
}
XV <- cbind(X,V)
for(h in 1:g)
{
	# A
	if(ncol(W) > 1)
	  A <- diag(1/c(W%*%c(sigma.e[,h])))
	else
	  A <- diag(1/c(W[,1]*sigma.e[1,h]))
	# B	
	B <- solve(sigma.b[,,h])
	# C
	C <- (1/sigma.c[h]) * diag(ncol(V))
	# E
	E <- solve(t(U) %*% A %*% U + B)
	#XVAU
	XVAU <- t(XV)%*%(A%*%U)
	# P
	P <- t(XV) %*% A %*% XV
	P[2+1:m,2+1:m] <- P[2+1:m,2+1:m] + C
	P <- P - XVAU %*% E %*% t(XVAU)
	P <- solve(P)/nh[h]
	#KOK
	XVAUEK <- XVAU %*% E %*% K2
	KOK <- (t(K1)-t(XVAUEK)) %*% P %*% K1
	KOK <- KOK - t(K1)%*% P %*%  XVAUEK + K2 %*% E %*% K2
	KOK <- KOK + t(XVAUEK)  %*% P %*%  XVAUEK
	omega[h] <- c(KOK)
}
sqrt(omega)	
}

scores.wire <-function(obj,contrast=NULL) 
{
if(obj$nb !=2 && obj$nb !=3) 
stop("only two or three classes can be compared.")

# get the sqrt KOK 
ooo <- eq8.wire(obj$m,obj$g,obj$nb, obj$X, obj$W, obj$U, obj$V, obj$sigma.e, obj$sigma.b, obj$sigma.c,colSums( obj$tau),contrast)
# contrasts
if(obj$nb ==2)
{
  if(is.null(contrast))
    d1d2 = (t(( obj$eb[,1,]- obj$eb[,2,]))+ ( obj$beta[1,]- obj$beta[2,]))/ooo
  else 
    d1d2 = (t(( obj$eb[,1,]  * contrast[1] + obj$eb[,2,]  * contrast[2] ))
           +  ( obj$beta[1,] * contrast[1] + obj$beta[2,] * contrast[2] ))/ooo  
} else { 
  
  if(is.null(contrast))
    d1d2 = (t(( obj$eb[,1,]- obj$eb[,3,]))+ ( obj$beta[1,]- obj$beta[3,]))/ooo
  else 
    d1d2 = (t(( obj$eb[,1,] * contrast[1] + obj$eb[,2,]  * contrast[2]  + obj$eb[,3,]  * contrast[3] ))
          +  ( obj$beta[1,] * contrast[1] + obj$beta[2,] * contrast[2]  + obj$beta[3,] * contrast[3] ))/ooo
}
# equation 5
c(rowSums(  obj$tau * t(d1d2)))
}


# permutation and null distribution

# B=99 permutations for class labels
# when calculate the W_j, only re-do the numerator,
# but keep the denominator same!

wj2.permuted <- function(data,ret,nB=99,contrast=NULL) {

mat.ABC.wire<-function(U,VV,W,sigma.e,DU,sigma.c,g,m)
{
		A<-B<-C<-BC<-array(0,dim=c(m,m,g))
		for(h in 1:g)
		{
			if(ncol(W)>1)
				A[,,h]<-diag(as.vector(W%*%sigma.e[,h]))
			else
			{
				A[,,h]<-diag(c(W[,1]*sigma.e[1,h]))
				
			}
			B[,,h]<-A[,,h]+U%*%DU[,,h]%*%t(U)
			C[,,h]<-sigma.c[h]*VV
			BC[,,h]<-B[,,h]+C[,,h]
		}
		list(A=A,B=B,C=C,BC=BC)
}
matM.wire<-function(B,C,tau,g,m)
{
		M<-array(0,dim=c(m,m,g))
		for(h in 1:g)
			M[,,h]<-solve(B[,,h]+C[,,h]*sum(tau[,h]))
		M
}
#--------------------------------------------
	
eb.estep<-function(tau,y,mu,DU,U,B,C,M,g,n,m,qb)
{
		eb1<-array(0,dim=c(n,qb,g))
		eb2<-array(0,c(g,qb,qb))
		invB<-array(0,c(m,m,g))
		
		for(h in 1:g)
		{
			invB[,,h]<-solve(B[,,h])
			# E(bi|y)
			
			eb1[,,h]<-t( DU[,,h]%*%t(U)%*%invB[,,h]%*%(  (t(y)-mu[,h])
																										-c(C[,,h]%*%M[,,h]%*%colSums(t(t(y)-mu[,h])*tau[,h])) ))
			#E(bi%*%bi^t|y)
			eb2[h,,]<-( sum(tau[,h])*DU[,,h]-(sum(tau[,h])-1)*DU[,,h]%*%t(U)%*%invB[,,h]%*%U%*%DU[,,h]	
									-DU[,,h]%*%t(U)%*%M[,,h]%*%U%*%DU[,,h])
	DU[,,h]<-(eb2[h,,]+ t(eb1[,,h]*tau[,h])%*%eb1[,,h])/sum(tau[,h])
}
		
list(DU=DU,eb1=eb1,invB=invB)
}
g  <- ret$g
X<-as.matrix(ret$X)	
U<-as.matrix(ret$U)
V<-as.matrix(ret$V)
W<-as.matrix(ret$W)
	
qb<-ncol(U)
qc<-ncol(V)
qe<-ncol(W)
	
VV<-V%*%t(V)
	
data<-as.matrix(data)

n<-nrow(data)
m<-ncol(data)
mu<-matrix(0,ncol=g,nrow=m)
	
for(h in 1:g) 
 mu[,h] <- as.vector(X %*% ret$beta[,h])

mobj <- mat.ABC.wire(U,VV,W,ret$sigma.e,ret$sigma.b,ret$sigma.c,g,m)
A <- mobj$A
B <- mobj$B
C <- mobj$C
BC<- mobj$BC	
#-------------------------------------
# get the sqrt KOK 

ooo <- eq8.wire(m,g,ret$nb,X,W,U,V,ret$sigma.e,ret$sigma.b, ret$sigma.c,colSums(ret$tau),contrast)
	
ooo[is.na(ooo) ] <- 1e+10
	
#-------------------------------------	
set.seed(1234)
wj0 <- array(0,c(n,nB))
for(b in 1:nB) { #do B permutation
da <- data[,sample(1:m,m,replace=FALSE)]	
#   get tau for da		
tau  <- tau.estep.wire(da,ret$pro,mu,BC,n,m,g)$tau
M <- matM.wire(B,C,tau,g,m)
###########################################
# update eb,			
eb <- eb.estep(tau,da,mu,ret$sigma.b,U,B,C,M,g,n,m,qb)$eb1
# ec is fixed
ec <- ret$ec
		
###########################################
		
# get new tau		
mu2    <- mu	
for(h in 1:g)	
mu2[,h]<- c(X%*%ret$beta[,h]+V%*%ec[,h])
tau    <- tau.estep.wire(da,ret$pro,mu2,B,n,m,g)$tau		
		
###########################################
				
# contrasts
if(ret$nb ==2)
{
  if(is.null(contrast))
    d1d2 = (t(( eb[,1,]- eb[,2,]))+ ( ret$beta[1,]- ret$beta[2,]))/ooo
  else 
    d1d2 = (t(( eb[,1,]  * contrast[1] + eb[,2,]  * contrast[2] ))
           +  ( ret$beta[1,] * contrast[1] + ret$beta[2,] * contrast[2] ))/ooo  
} else { 
  
  if(is.null(contrast))
    d1d2 = (t(( eb[,1,]- eb[,3,]))+ ( ret$beta[1,]- ret$beta[3,]))/ooo
  else 
    d1d2 = (t(( eb[,1,] * contrast[1] + eb[,2,]  * contrast[2]  + eb[,3,]  * contrast[3] ))
          +  ( ret$beta[1,] * contrast[1] + ret$beta[2,] * contrast[2]  + ret$beta[3,] * contrast[3] ))/ooo
}
# equation 5
wj0[,b] <- (rowSums(  tau * t(d1d2)))
		
} #end of B permutations	
wj0
}

pvalue.wire <- function(wj,wj0){	
n0 <- length(wj)
nn <- length(wj0)
pv <- rep(0,n0)
for(j in 1:n0) {
pv[j] <- sum( abs(c(wj,wj0)) >= abs(wj[j])) /(nn+n0)
}
pv	
}

#---------------------------------
# end of program
#---------------------------------

#all normal distRand versions now moved to C and to the test files
#tdist <-function(nu,U,z.list,mu,gamm){
#	#use nu and z.list to get D which is scale matrix
#	eek<-getEk(z.list)
#	Dvecs<-Map("*",eek,nu)
#	Dvec<-addVecs(Dvecs) #at this point still a vector
#	Dmat<-diag(Dvec)

#	value<-dmvt(U,delta=mu,sigma=Dmat,log=TRUE,df=gamm,type="shifted")
#	list(value=value)

#}

#note that Dstarinv is just the DIAGONAL of Dstarinv
tdist2<-function(tconst,u, Dstarinv,zeta,myq){
	inside<-1+sum(t(u)*Dstarinv*u)/zeta
	logft<-tconst - ((zeta+myq)/2)*log(inside)
	logft
}

tconstant<-function(zeta,myq,Dstarinvdiag){
	piece1<-log(gamma((zeta+myq)/2))
	piece1b<--log(gamma(zeta/2))
	piece1c<--myq*log(zeta*pi)/2
	piece2<- .5*log(prod(Dstarinvdiag))	
	piece1+piece1b+piece1c+piece2
}

# no longer used
#distRandGeneral2<-function(uvec,mu,Sigma.inv,logDetSigmaInv){
#	print("still using distRandGeneral2.R")
#	umu<-uvec-mu
#	piece2<-t(umu)%*%Sigma.inv%*%umu

#	as.vector(.5*(logDetSigmaInv-piece2))
#}

# old version of distRand, used only for tests





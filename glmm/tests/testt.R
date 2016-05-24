# test the t distribution

library(glmm)

tdist2<-function(tconst,u, Dstarinv,zeta,myq){
	inside<-1+t(u)%*%Dstarinv%*%u/zeta
	logft<-tconst - ((zeta+myq)/2)*log(inside)
	as.vector(logft)
}

zeta<-5
myq<-10
set.seed(1234)
u<-rt(10,df=zeta)
Dstarinvdiag<-rep(1.1,10)
tconstant<-glmm:::tconstant
tconst<-tconstant(zeta,myq,Dstarinvdiag)

Dstarinv<-diag(Dstarinvdiag)
tR<-tdist2(tconst,u,Dstarinv,zeta,myq)

#compare version in R to version in C
tC<-.C("tdist",as.double(Dstarinv),  as.integer(myq), as.double(u), as.integer(zeta), as.double(tconst), double(1))
all.equal(tC[[6]],tR)

#compare my version to the mvtnorm version
Dstar<-diag(10)/1.1
theirs<-dmvt(u,sigma=Dstar,df=zeta)
theirs==tR

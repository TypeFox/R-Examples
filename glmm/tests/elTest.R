#going to test el for the Booth and Hobert example

library(glmm)
data(BoothHobert)
set.seed(1234)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobert, family.glmm=bernoulli.glmm, m=21, doPQL=TRUE, debug=TRUE)

mod.mcml<-mod.mcml1$mod.mcml
z<-mod.mcml$z[[1]]
x<-mod.mcml$x
y<-mod.mcml$y

stuff<-mod.mcml1$debug
beta.pql<-stuff$beta.pql
nu.pql<-stuff$nu.pql
u.pql<-u.star<-stuff$u.star
umat<-stuff$umat

family.glmm<-bernoulli.glmm

objfun<-glmm:::objfun
getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

############################################
#this should be the same as el
getFamily<-glmm:::getFamily
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

#compare elR and el.C for a value of eta
eta<-rep(2,150)
that<-elR(mod.mcml$y,mod.mcml$x,eta,family.mcml=bernoulli.glmm)	
this<-.C("elc",as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(1),value=double(1),gradient=double(ncol(mod.mcml$x)),hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.numeric(that$value),this$value)
all.equal(as.numeric(that$gradient),this$gradient)
all.equal(as.numeric(that$hessian),this$hessian)

#compare to elval
elvalout<-.C("elval",as.double(mod.mcml$y),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(1),value=double(1))
all.equal(as.numeric(that$value),elvalout$value)

#compare to elGH
elGHout<-.C("elGH",as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(1),gradient=double(ncol(mod.mcml$x)),hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.numeric(that$gradient),elGHout$gradient)
all.equal(as.numeric(that$hessian),elGHout$hessian)




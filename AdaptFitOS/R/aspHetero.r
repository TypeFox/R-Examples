aspHetero=function(object,xx,nknots=5,knots=NULL,basis="os", degree=c(3,2),tol=1e-8,niter=100,niter.var=250){  #  ,basis="os",degree=3
data.fit=object

y=object$info$y

fit.save=data.fit
fit.start =fit.save
data.fit=fit.start
#
fit<-fitted(data.fit)
#

n=length(object$fitted)
all<-rep(1,n)
X<-data.fit$design.matrices$Xb
Z<-data.fit$design.matrices$Zb
C<-data.fit$design.matrices$Wb
K<-ncol(Z)
D<-diag(c(rep(0,ncol(X)),rep(1,K)))

sigma.eps<-data.fit$y.cov[1]
sigma.u<-data.fit$random.var
beta<-data.fit$fit$coef$fixed


nk=nknots
if (basis=="os"){
  knotsx= seq(min(xx),max(xx),length=nk+2)[-c(1,nk+2)]
  smooth=list(m=degree,bs.dim=nk,term="xx",p.order=degree,knots=knotsx)
  class(smooth)= "ospline.smooth"
} else {
if (is.null(knots)) kn.e= seq(min(xx),max(xx),length=nk)
else {kn.e=knots; kn=length(knots)}
names(kn.e) <- NULL
  smooth=list(m=c(degree,degree+1),bs.dim=nk,knots= kn.e,term="xx")
  if (basis=="trunc.poly") class(smooth)="tlspline.smooth"
  else class(smooth)= basis
}
data.xx = data.frame(xx = xx)
CZ.temp <- Predict.matrix.lme(smooth,data.xx,center=T)
if (basis=="os") Xe <- CZ.temp$C          
else Xe <- CZ.temp$C[,-1,drop=F]
Ze <- CZ.temp$Z

Ce=cbind(Xe,Ze)


fit.start$sigmax=list(p.order=degree,m=degree,bs.dim=nk,knots= knotsx,term="xx")
class(fit.start$sigmax)=  data.fit$info$pen$basis

De<-(c(rep(0,ncol(Xe)),rep(1,ncol(Ze))))
var.fit<-log((((y-fit)-mean(y-fit))^2)/var(y-fit))
omega.fit<-lme(var.fit~Xe-1,random=list(all=pdIdent(~Ze-1)))
omega<-c(omega.fit$coef$fixed,as.vector(unlist(omega.fit$coef$random)))
sigma.v<- as.vector(omega.fit$sigma^2 * exp(2 * unlist(omega.fit$modelStruct$reStruct)))
sigma.v<-sqrt(sigma.v)


m<-1
for(j in 1:niter){
    res<-c(y-X%*%beta)
    l<-1
    for(i in 1:niter.var) {
        Sigma.eps<-exp(c(Ce%*%omega))
        Sigma.eps.inv<-exp(-c(Ce%*%omega))
        Zi <- ginv(t(Z)%*%(Z*Sigma.eps.inv)/sigma.eps+diag(K)/sigma.u)%*%t(Z*Sigma.eps.inv)/sigma.eps
        A <- apply(Z*(t(Zi)),1,sum)
        AA <- rep(1,n)-2*(A)+apply(Z*(t((Zi%*%Z)%*%Zi)),1,sum)
        e <- c((Sigma.eps.inv)*res/sigma.eps-(Sigma.eps.inv)*(Z%*%(Zi%*%res)/sigma.eps))

        cac<-t(Ce)%*%((AA/2)*Ce)
        u<- cac%*%omega+t(Ce)%*%(Sigma.eps*(e^2)*sigma.eps-rep(1,n)+(A))/2
        Ivv<-ginv(cac+diag(De)/sigma.v)
        omega1<- Ivv%*%u

        epsilon.omega<-sum((omega-omega1)^2)/sum(omega^2)
        l<-l+1
        if (epsilon.omega <= tol)  break
        if (l == niter.var) stop("Iteration limit reached without convergence in variance parameter")
        omega<-omega1
    }
    zaz<-t(Ze)%*%((AA/2)*Ze)
    IvvZ<-solve(zaz+diag(ncol(zaz))/sigma.v)
    sigma.v<-c(t(omega)%*%(De*omega)/sum(diag(zaz%*%IvvZ)))

    y1<-(exp(-c(Ce%*%omega)/2))*y
    X1<-(exp(-c(Ce%*%omega)/2))*X
    Z1<-(exp(-c(Ce%*%omega)/2))*Z

    data.fit<-lme(y1~X1-1,random=list(all=pdIdent(~Z1-1)))
 #
  fit1<- (exp(c(Ce%*%omega)/2))*fitted(data.fit)
 #
   sigma.eps<-data.fit$sigma^2
    sigma.u<-as.vector(data.fit$sigma^2 * exp(2 * unlist(data.fit$modelStruct$reStruct)))
    beta<-as.vector(data.fit$coef$fixed)
    uhat= c(data.fit$coef$random$all)
    epsilon.fit <- sum((fit -fit1)^2)/sum(fit^2)
    fit <- fit1
    m<-m+1

    if (epsilon.fit <= tol) break
    if (m == niter) stop("Iteration limit reached without convergence")
}

fit.start$fitted=fit
fit.start$coef.mean=c(beta,uhat)
fit.start$fit$sigma         =  sqrt( sigma.eps)

fit.start$sigmax$fitted=sqrt(Sigma.eps)
fit.start$sigmax$coef=omega
return(fit.start)
}

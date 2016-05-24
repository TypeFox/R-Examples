estimateParameters.copula <- function(object,...) {

trend<-NULL
if(is.null(object$formulaString)){
  F<-as.matrix(rep(1,length(object$observations@data$value)))
  object$formulaString<-value~1
} else {
  F<-model.matrix(object$formulaString,data=object$observations)
}
trend$F<-F
depVar = as.character(object$formulaString[[2]])


# If subsampling of data is necessary set subsample<Inf
ndim = dim(object$observations@data)[1]
subsample<-Inf
if(ndim >subsample){
  index<-sample(1:ndim,subsample)
} else {
  index<-1:ndim
}

newdata = object$observations[index,]
trend$F<-as.matrix(trend$F[index,])

anisotropy<-NULL
if(!"anisPar" %in% names(object)){
object<-estimateAnisotropy(object)
}

anisotropy<-NULL
correlation<-NULL
if(object$anisPar$doRotation==T){
anisotropy$lower<-c(0,1)
anisotropy$upper<-c(pi,Inf)
if(is.null(object$copulaParams$anisotropy)){
anisotropy$params<-c(object$anisPar$direction*pi/180,object$anisPar$ratio)
} else {
anisotropy$params<-object$copulaParams$anisotropy$params
}
if(is.null(object$copulaParams$correlation)){
temp2_obj<-rotateAnisotropicData(object)
correlationmodel<-autofitVariogram(object$formulaString,temp2_obj,model=c("Ste","Sph","Gau","Exp"),kappa=c(0.2,1,2,5))
rm(temp2_obj)
} else {
correlation<-object$copulaParams$correlation
}
} else  {
if(is.null(object$copulaParams$correlation)){
correlationmodel<-autofitVariogram(object$formulaString,object$observations,model=c("Ste","Sph","Gau","Exp"),kappa=c(0.5,1,2,5))
} else {
correlation<-object$copulaParams$correlation
}
}

if(is.null(object$copulaParams$correlation)){
correlation$model<-as.character(correlationmodel$var_model$model[2])
correlation$params[1]<-max(0.01,correlationmodel$var_model$psill[1]/sum(correlationmodel$var_model$psill))
correlation$params[2]<-max(0.01,correlationmodel$var_model$range[2])
if(correlation$model=="Ste"){
correlation$params[3]<-correlationmodel$var_model$kappa[2]
}
}

if(correlation$model=="Ste"){
correlation$lower<-c(0.01,0.01,0.01)
correlation$upper<-c(0.99,Inf,20)
} else {
correlation$lower<-c(0.01,0.01)
correlation$upper<-c(0.99,Inf)
}
temp_obj<-object
temp_obj$predictionLocations<-NULL
temp_obj$observations<-newdata

margin<-NULL
if(is.null(object$copulaParams$margin$name)){
  marg<-findmargins(newdata@data[[depVar]])
if(length(marg$params)==2){
  margin$params<-marg$params[2]
} else {
  margin$params<-marg$params[2:3]
}
margin$name<-marg$name
trend$params<-marg$params[1]
marg<-findmarginbounds(margin)
if(length(margin$params)==1){
margin$lower<-marg$bounds$lower[2]
margin$upper<-marg$bounds$upper[2]
} else {
margin$lower<-marg$bounds$lower[2:3]
margin$upper<-marg$bounds$upper[2:3]
}
trend$lower<-marg$bounds$lower[1]
trend$upper<-marg$bounds$upper[1]
} else {
if(is.null(object$copulaParams$margin$params) || is.null(object$copulaParams$trend$params)){
marg<-findmargins(newdata@data[[depVar]],name=object$copulaParams$margin$name)
if(length(marg$params)==2){
margin$params<-marg$params[2]
} else {
margin$params<-marg$params[2:3]
}
trend$params<-marg$params[1]
margin$name<-marg$name
marg<-findmarginbounds(margin)
if(length(margin$params)==1){
margin$lower<-marg$bounds$lower[2]
margin$upper<-marg$bounds$upper[2]
} else {
margin$lower<-marg$bounds$lower[2:3]
margin$upper<-marg$bounds$upper[2:3]
}
trend$lower<-marg$bounds$lower[1]
trend$upper<-marg$bounds$upper[1]
} else {
if(is.null(object$copulaParams$margin$bounds)){
margin$name<-object$copulaParams$margin$name
trend$params<-object$copulaParams$trend$params
margin$params<-object$copulaParams$margin$params
marg<-findmarginbounds(list(name=object$copulaParams$margin$name,params=object$copulaParams$margin$params))
if(length(margin$params)==1){
margin$lower<-marg$bounds$lower[2]
margin$upper<-marg$bounds$upper[2]
} else {
margin$lower<-marg$bounds$lower[2:3]
margin$upper<-marg$bounds$upper[2:3]
}
trend$lower<-marg$bounds$lower[1]
trend$upper<-marg$bounds$upper[1]
}
}
}

if(dim(trend$F)[2]>1){
trend$lower<-c(trend$lower,rep(-Inf,dim(trend$F)[2]-1))
trend$upper<-c(trend$upper,rep(Inf,dim(trend$F)[2]-1))
if(is.null(object$copulaParams$trend$params)){
trend$params<-c(trend$params,rep(0,dim(trend$F)[2]-1))
}
}


copula<-NULL
if(is.null(object$copulaParams$copula$method)){
	copula$method="norm"
} else {
	copula<-object$copulaParams$copula
}
#copula$params<-2
#copula$lower<-0
#copula$upper<-Inf

estimation<-copulaEstimation(temp_obj,margin,trend,correlation,anisotropy,copula,...)
estimation$trend$F<-F
object$copulaParams<-estimation
object
}

copulaEstimation <- function(obj,margin,trend,correlation,anisotropy,copula,tol=0.001,...){
data<-NULL
data$x<-coordinates(obj$observations)[,1]
data$y<-coordinates(obj$observations)[,2]
depVar = as.character(obj$formulaString[[2]])
data$z<-obj$observations@data[[depVar]]

if (obj$params$debug.level > 1) cat(margin$name, "-Distribution", "\n")
distribfunction<-get(paste("p",margin$name,sep=""),mode="function")
densityfunction<-get(paste("d",margin$name,sep=""),mode="function")

if(copula$method=="chisq"){
params<-profilelikelihood(margin,correlation,trend,anisotropy,copula,data,distribfunction,densityfunction,tol,obj$params$debug.level,...)
} else {
if(copula$method=="norm"){
     params<-profilelikelihood2(margin,correlation,trend,anisotropy,data,distribfunction,densityfunction,tol,obj$params$debug.level,...)
}
}
params$copula$method<-copula$method
params
}

`findmargins` <-
function(data,name=NULL){
if(!is.null(name)){
margins<-name
} else {
if(sum(data<=0)>0){
margins<-c("norm","gev","t","logis")
} else {
margins<-c("norm","lnorm","gev","t","logis")
#add "gamma"
}
}
params<-NULL
name<-NULL
pvalue<-NULL
for(i in 1:length(margins)){
funct<-get(paste(margins[i],"testing",sep=""),mode="function")
estim<-funct(data)
if(i==1 || estim$pvalue<pvalue){#> wenn nach pvalue gefragt ist und < wenn nach statistic gefragt ist
name<-margins[i]
params<-estim$params
pvalue<-estim$pvalue
}
}
list(name=name,params=params)
}

`findmarginbounds` <-
function(margin){
switch(margin$name,
norm=list(name=margin$name,bounds=list(lower=c(-Inf,0.01),upper=c(Inf,Inf)),params=margin$params),
lnorm=list(name=margin$name,bounds=list(lower=c(-Inf,0.01),upper=c(Inf,Inf)),params=margin$params),
t=list(name=margin$name,bounds=list(lower=c(0.01,-Inf),upper=c(Inf,Inf)),params=margin$params),
gev=list(name=margin$name,bounds=list(lower=c(-Inf,0.01,-Inf),upper=c(Inf,Inf,Inf)),params=margin$params),
gamma=list(name=margin$name,bounds=list(lower=c(0.01,0.01),upper=c(Inf,Inf)),params=margin$params),
logis=list(name=margin$name,bounds=list(lower=c(-Inf,0.01),upper=c(Inf,Inf)),params=margin$params)
)
}

`gammatesting` <-
function(data){
params<-try(fitdistr(data,"gamma")$estimate,TRUE)
if(class(params)!="try-error"){
pvalue<-ks.test(data,pgamma,params[1],params[2])$statistic
l<-list(params=params,pvalue=pvalue)
} else {
l<-list(params=NULL,pvalue=Inf)
}
l
}

`logistesting` <-
function(data){
params<-try(fitdistr(data,"logistic")$estimate,TRUE)
if(class(params)!="try-error"){
pvalue<-ks.test(data,plogis,params[1],params[2])$statistic
l<-list(params=params,pvalue=pvalue)
} else {
l<-list(params=NULL,pvalue=Inf)
}
l
}

`gevtesting` <-
function(data){
#params<-try(fgev(data,std.err=FALSE)$param,TRUE)
params<-gevfit(data)
#if(class(params)!="try-error"){
if(sum(is.nan(params))==0){
	pvalue<-ks.test(data,pgev,params[1],params[2],params[3])$statistic
	l<-list(params=params,pvalue=pvalue)
} else {
	l<-list(params=NULL,pvalue=Inf)
}
l
}

`lnormtesting` <-
function(data){
normtesting(log(data))
}

`ttesting` <-
function(data){
params<-try(fitdistr(data,"t",start=NULL,s=1)$estimate,TRUE)
if(class(params)!="try-error"){
pvalue<-ks.test(data,pt,params[2],params[1])$statistic
l<-list(params=c(params[2],params[1]),pvalue=pvalue)
} else {
l<-list(params=NULL,pvalue=Inf)
}
l
}


`normtesting` <-
function(data){
mu<-sum(data)/length(data)
sigma<-sqrt(1/(length(data)-1)*sum((data-mu)^2))
pvalue<-ks.test(data,pnorm,mu,sigma)$statistic
list(params=c(mu,sigma),pvalue=pvalue)
}


`likelianiso` <-
function(margin,trend,correlation,anisotropy,data,h,distr,dens,copula,anisomode,debug.level){

if(anisomode){
      xy<-matrix(c(cos(anisotropy$params[1]),sin(anisotropy$params[1]),-anisotropy$params[2]*sin(anisotropy$params[1]),anisotropy$params[2]*cos(anisotropy$params[1])),ncol=2,byrow=TRUE)%*%rbind(t(data$x),t(data$y))
	x<-xy[1,]
	y<-xy[2,]
	h<-as.matrix(dist(t(as.matrix(rbind(x,y)))))
}

expected<-trend$F%*%trend$params;
if(length(margin$params)==1){
z<-distr(data$z,expected,margin$params[1])
zz<-dens(data$z,expected,margin$params[1],log=TRUE)
} else {
z<-distr(data$z,expected,margin$params[1],margin$params[2])
zz<-dens(data$z,expected,margin$params[1],margin$params[2],log=TRUE)
}
ok<-(sum(z!=1)==length(z))*(sum(z!=0)==length(z));
lambda<-copula$params
if(lambda>=0 & ok){
lik<-0
z<-qchisq(z,1,lambda)
x<-dchisq(z,1,lambda)
z<-sqrt(z)
lambda<-sqrt(lambda)
m<-c(lambda,lambda)
cov=covar(h,correlation$params,correlation$model)
for(i in 1:(length(data$z)-1)){
    for(j in (i+1):length(data$z)){            
		lik<-lik-(log(sum(dmvnorm(matrix(c(z[i],z[j],-z[i],z[j],z[i],-z[j],-z[i],-z[j]),ncol=2,byrow=TRUE),m,cov[c(i,j),c(i,j)])))-log(z[i])-log(z[j])-log(x[i])-log(x[j]))-zz[i]-zz[j]
	}
}
} else {
lik<-10^10
}
if(debug.level>1) cat(lik, c(margin$params,trend$params,correlation$params,anisotropy$params,copula$params),"\n")
lik
}

`likelianiso2` <-
function(margin,trend,correlation,anisotropy,data,h,distr,dens,anisomode,sigmainv,logSqrtDetSigma,debug.level){

if(anisomode){
      xy<-matrix(c(cos(anisotropy$params[1]),sin(anisotropy$params[1]),-anisotropy$params[2]*sin(anisotropy$params[1]),anisotropy$params[2]*cos(anisotropy$params[1])),ncol=2,byrow=TRUE)%*%rbind(t(data$x),t(data$y))
	x<-xy[1,]
	y<-xy[2,]
	h<-as.matrix(dist(t(as.matrix(rbind(x,y)))))
}

expected<-trend$F%*%trend$params;
if(length(margin$params)==1){
	z<-distr(data$z,expected,margin$params[1])
	zz<-dens(data$z,expected,margin$params[1],log=TRUE)
} else {
	z<-distr(data$z,expected,margin$params[1],margin$params[2])
	zz<-dens(data$z,expected,margin$params[1],margin$params[2],log=TRUE)
}
ok<-(sum(z!=1)==length(z))*(sum(z!=0)==length(z))
if(ok==1){
	z<-qnorm(z,0,1)
	if(is.null(sigmainv) || is.null(logSqrtDetSigma)){
		cov<-covar(h,correlation$params,correlation$model)
		#sing<-svd(cov)
		#detK<-sum(log(sing[[1]]))
		#invK<-sing[[2]]%*%diag(1/sing[[1]])%*%t(sing[[2]])
		#lik<-0.5*detK+ 0.5*t(z)%*%invK%*%z-0.5*t(z)%*%z-sum(zz)
		lik<--logdmvnormnoconst(z,cov)-0.5*t(z)%*%z-sum(zz)
	} else {
		lik<-0.5*t(z)%*%sigmainv%*%z+logSqrtDetSigma-0.5*t(z)%*%z-sum(zz)
	}
} else {
lik<-10^5
}
if(debug.level>1) cat(lik, c(margin$params,trend$params,correlation$params,anisotropy$params),"\n")
lik
}

`profilelikelihood` <-
function(margin,correlation,trend,anisotropy,copula,data,distr,dens,tolerance,debug.level,...){

if(is.null(anisotropy$params)){
h<-as.matrix(dist(t(as.matrix(rbind(data$x,data$y)))))
} else {
        xy<-matrix(c(cos(anisotropy$params[1]),sin(anisotropy$params[1]),-anisotropy$params[2]*sin(anisotropy$params[1]),anisotropy$params[2]*cos(anisotropy$params[1])),ncol=2,byrow=TRUE)%*%rbind(t(data$x),t(data$y))
x<-xy[1,]
y<-xy[2,]
h<-as.matrix(dist(t(as.matrix(rbind(x,y)))))
    }


    len=length(margin$params);
    lik=1000000;
    oldlik=0;


    while(abs(lik-oldlik)>tolerance){
        oldlik=lik;
        if(!is.null(correlation$params)){
            res<-optim(correlation$params,optimfun5,gr=NULL,margin=margin,trend=trend,correlation=correlation,anisotropy=anisotropy,data=data,h=h,distr=distr,dens=dens,copula=copula,anisomode=0,debug.level,method="L-BFGS-B",lower=correlation$lower,upper=correlation$upper,...)
correlation$params<-res$par
lik<-res$value
        }
       if(!is.null(margin$params) || !is.null(trend$params)){
res<-optim(c(margin$params,trend$params),optimfun4,gr=NULL,margin=margin,trend=trend,correlation=correlation,anisotropy=anisotropy,data=data,h=h,len=len,distr=distr,dens=dens,copula=copula,anisomode=0,debug.level,method="L-BFGS-B",lower=c(margin$lower,trend$lower),upper=c(margin$upper,trend$upper),...)
margin$params<-res$par[1:len]
trend$params<-res$par[(len+1):length(res$par)]
lik<-res$value
        }

if(!is.null(copula$params)){
            res<-optim(copula$params,optimfun7,gr=NULL,margin=margin,trend=trend,correlation=correlation,anisotropy=anisotropy,data=data,h=h,distr=distr,dens=dens,copula=copula,anisomode=0,debug.level,method="L-BFGS-B",lower=copula$lower,upper=copula$upper,...)
copula$params<-res$par
lik<-res$value
        }

        if(!is.null(anisotropy$params)){
            res<-optim(anisotropy$params,optimfun6,gr=NULL,margin=margin,trend=trend,correlation=correlation,anisotropy=anisotropy,data=data,h=h,distr=distr,dens=dens,copula=copula,anisomode=1,debug.level,method="L-BFGS-B",lower=anisotropy$lower,upper=anisotropy$upper,...)
anisotropy$params<-res$par
lik<-res$value
        xy<-matrix(c(cos(anisotropy$params[1]),sin(anisotropy$params[1]),-anisotropy$params[2]*sin(anisotropy$params[1]),anisotropy$params[2]*cos(anisotropy$params[1])),ncol=2,byrow=TRUE)%*%rbind(t(data$x),t(data$y))
x<-xy[1,]
y<-xy[2,]
h<-as.matrix(dist(t(as.matrix(rbind(x,y)))))
}
    }
    optimized<-NULL
    optimized$margin=margin
    optimized$correlation=correlation
    optimized$anisotropy=anisotropy
    optimized$trend=trend
optimized$copula=copula
optimized
}

`profilelikelihood2` <-
function(margin,correlation,trend,anisotropy,data,distr,dens,tolerance,debug.level,...){

if(is.null(anisotropy$params)){
	h<-as.matrix(dist(t(as.matrix(rbind(data$x,data$y)))))
} else {
	h<-NULL
}


    len=length(margin$params);
    lik=1000000;
    oldlik=0;

    while(abs(lik-oldlik)>tolerance){
        oldlik=lik;
	  if(!is.null(anisotropy$params)){
        	xy<-matrix(c(cos(anisotropy$params[1]),sin(anisotropy$params[1]),-anisotropy$params[2]*sin(anisotropy$params[1]),anisotropy$params[2]*cos(anisotropy$params[1])),ncol=2,byrow=TRUE)%*%rbind(t(data$x),t(data$y))
	  	x<-xy[1,]
	  	y<-xy[2,]
	  	h<-as.matrix(dist(t(as.matrix(rbind(x,y)))))
	  }
	  cov=covar(h,correlation$params,correlation$model)
        R=chol(cov)
    	  sigmainv <- chol2inv(R)
    	  logSqrtDetSigma <- sum(log(diag(R)))
	  

       if(!is.null(margin$params) || !is.null(trend$params)){
		res<-optim(c(margin$params,trend$params),optimfun1,gr=NULL,margin=margin,trend=trend,correlation=correlation,anisotropy=anisotropy,data=data,h=h,len=len,distr=distr,dens=dens,anisomode=0,sigmainv,logSqrtDetSigma,debug.level,method="L-BFGS-B",lower=c(margin$lower,trend$lower),upper=c(margin$upper,trend$upper),...)
		margin$params<-res$par[1:len]
		trend$params<-res$par[(len+1):length(res$par)]
		lik<-res$value
        }

        if(!is.null(correlation$params)){
            res<-optim(correlation$params,optimfun2,gr=NULL,margin=margin,trend=trend,correlation=correlation,anisotropy=anisotropy,data=data,h=h,distr=distr,dens=dens,anisomode=0,debug.level,method="L-BFGS-B",lower=correlation$lower,upper=correlation$upper,...)
		correlation$params<-res$par
		lik<-res$value
        }

        if(!is.null(anisotropy$params)){
            res<-optim(anisotropy$params,optimfun3,gr=NULL,margin=margin,trend=trend,correlation=correlation,anisotropy=anisotropy,data=data,h=h,distr=distr,dens=dens,anisomode=1,debug.level,method="L-BFGS-B",lower=anisotropy$lower,upper=anisotropy$upper,...)
		anisotropy$params<-res$par
		lik<-res$value
        	#xy<-matrix(c(cos(anisotropy$params[1]),sin(anisotropy$params[1]),-anisotropy$params[2]*sin(anisotropy$params[1]),anisotropy$params[2]*cos(anisotropy$params[1])),ncol=2,byrow=TRUE)%*%rbind(t(data$x),t(data$y))
		#x<-xy[1,]
		#y<-xy[2,]
		#h<-as.matrix(dist(t(as.matrix(rbind(x,y)))))
	  }
    }
    optimized<-NULL
    optimized$margin=margin;
    optimized$correlation=correlation;
    optimized$anisotropy=anisotropy;
    optimized$trend=trend;
optimized
}

`optimfun1` <-
function(params,margin,trend,correlation,anisotropy,data,h,len,distr,dens,anisomode,sigmainv,logSqrtDetSigma,debug.level){
        margin$params=params[1:len];
        trend$params=params[(len+1):length(params)];
        likelianiso2(margin,trend,correlation,anisotropy,data,h,distr,dens,anisomode,sigmainv,logSqrtDetSigma,debug.level)
}

`optimfun2` <-
function(params,margin,trend,correlation,anisotropy,data,h,distr,dens,anisomode,debug.level){
        correlation$params=params;
        likelianiso2(margin,trend,correlation,anisotropy,data,h,distr,dens,anisomode,NULL,NULL,debug.level)
}

`optimfun3` <-
function(params,margin,trend,correlation,anisotropy,data,h,distr,dens,anisomode,debug.level){
        anisotropy$params=params;
        likelianiso2(margin,trend,correlation,anisotropy,data,h,distr,dens,anisomode,NULL,NULL,debug.level)
}

`optimfun4` <-
function(params,margin,trend,correlation,anisotropy,data,h,len,distr,dens,copula,anisomode,debug.level){
        margin$params=params[1:len];
        trend$params=params[(len+1):length(params)];
        likelianiso(margin,trend,correlation,anisotropy,data,h,distr,dens,copula,anisomode,debug.level)
}

`optimfun5` <-
function(params,margin,trend,correlation,anisotropy,data,h,distr,dens,copula,anisomode,debug.level){
        correlation$params=params;
        likelianiso(margin,trend,correlation,anisotropy,data,h,distr,dens,copula,anisomode,debug.level)
}

`optimfun6` <-
function(params,margin,trend,correlation,anisotropy,data,h,distr,dens,copula,anisomode,debug.level){
        anisotropy$params=params;
        likelianiso(margin,trend,correlation,anisotropy,data,h,distr,dens,copula,anisomode,debug.level)
}

`optimfun7` <-
function(params,margin,trend,correlation,anisotropy,data,h,distr,dens,copula,anisomode,debug.level){
        copula$params=params;
        likelianiso(margin,trend,correlation,anisotropy,data,h,distr,dens,copula,anisomode,debug.level)
}



`logdmvnormnoconst`<-
function (x, sigma) {
    R<-chol(sigma)
    sigmainv <- chol2inv(R)
    logdet <- sum(log(diag(R)))
    -logdet-0.5*t(x)%*%sigmainv%*%x 
}




`gevfit`<-function(x,...){

n = length(x);
x = sort(x);
xmin = x[1];
xmax = x[n];
rangex = xmax-xmin;

if(n == 0 || !is.finite(rangex)){
    parmhat = NaN;
    return(parmhat)
} else {
	if(rangex < 1e-10){
    parmhat = c(0, 0, x[1]);
    return(parmhat)
	}
}

F = (0.5:(n-0.5))/ n;

init<-function(k,F,x) 1-cor(x,qgev(F,0,1,k))
k0 = optim(0,init,F=F,x=x,method="BFGS")$par;

b<-lsfit(qgev(F,0,1,k0),x)$coefficients
sigma0 = b[2];
mu0 = b[1];
if((k0 < 0 && (xmax > -sigma0/k0+mu0)) || (k0 > 0 && (xmin < -sigma0/k0+mu0))){

    p<-try(fgev(data,std.err=FALSE)$param)
    if(class(p)!="try-error"){
    	 mu0=p[1]
    	 sigma0=p[2]
       k0=0;
    } else {
	 return(NaN)
    }
}
parmhat = c(k0,log(sigma0),mu0);

output = optim(parmhat,negloglike,data=x,control=list(reltol=1e-10,maxit=10000))$par

parmhat[2] = exp(output[2]);
parmhat[1] = output[3]
parmhat[3] = output[1]
parmhat
}


`negloglike`<-function(parms, data){

k     = parms[1];
lnsigma = parms[2];
sigma = exp(lnsigma);
mu    = parms[3];

n = length(data);
z = (data - mu) / sigma;

if(abs(k) > 1e-16){
    t = 1 + k*z
    if(min(t) > 0){
        u = 1 + k*z
        lnu = log1p(k*z); # log(1 + k*z)
        t = exp(-(1/k)*lnu); # (1 + k*z).^(-1/k)
        nll = n*lnsigma + sum(t) + (1+1/k)*sum(lnu)
    } else {
        nll = 1e10
    }
} else {
    nll = n*lnsigma + sum(exp(-z) + z)

}

nll
}



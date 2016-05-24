fitPP.fun <-
function(covariates=NULL,start, fixed = list(),  
posE=NULL,inddat=NULL,POTob=NULL,nobs=NULL,
tind=TRUE, tim=NULL, minfun="nlminb",modCI=TRUE, CIty="Transf", clevel=0.95,
tit="", modSim="FALSE",dplot=TRUE,xlegend="topleft",lambdaxlim=NULL,
lambdaylim=NULL, ...)
{

call<-match.call()
if ( is.null(posE)&is.null(POTob) )
{
stop('Error: At least one of the arguments, posE or POTob must be specified')
}

if (!is.list(fixed))     stop("'fixed' must be a named list")
if (!is.list(start) || is.null(names(start)) )     stop("'start' must be a named list")
if (any(!names(fixed) %in% names(start))) 
        stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
 

if (is.null(POTob)==FALSE) 
{
T<-POTob$T
thres<-POTob$thres
nobs<-length(T)


exc<-as.numeric(T>thres)
inrachtx<-(c(0,diff(exc))==1)
numerachtx<-cumsum(inrachtx)[exc==1]
intentx<-T[exc==1]
posicintx<-c(1:nobs)[inrachtx==1]
posE<-posicintx+tapply(intentx,INDEX=numerachtx, FUN=which.max)-1

inddat<-1-exc
inddat[posE]<-1
}

posE<-as.numeric(posE)

if (!is.null(covariates)) 
{	covariates<-as.matrix(covariates)
	if (min(dim(covariates))==0) covariates<-NULL
}


if ( is.null(nobs)&is.null(covariates)&is.null(inddat) )
{
stop('Error: At least one of the arguments, nobs, covariates,inddat or POTob must be specified')
}
if (is.null(nobs)& is.null(inddat)) nobs<-dim(covariates)[1]
if (is.null(nobs)& is.null(covariates)) nobs<-length(inddat)

if (is.null(inddat)) inddat<-rep(1,nobs)
nobs<-length(inddat)
control<-sum(inddat)

if ((tind==FALSE)&(is.null(covariates)==TRUE))
{
stop('Error: Model without covariates and without constant term')
}

if (!is.null(covariates))
{
	if (nobs!=dim(covariates)[1]) stop('Error:  Number of covariate observations  is not equal to nobs')
}

namcovariates<-dimnames(covariates)[[2]]
ncov<-dim(covariates)[2]
if (is.null(namcovariates)&(!is.null(ncov))) namcovariates<- paste('Covariate',c(1:ncov),sep='')
if (tind==TRUE)
{
	covariates<-cbind(rep(1,nobs),covariates)
      namcovariates<-c('Intercept',namcovariates)
}
dimnames(covariates)<-list(NULL,namcovariates)


if (is.null(tim)==TRUE) tim<-c(1:nobs)


Bcovariates<-as.matrix(covariates[inddat==1,])
covariatest<-as.matrix(covariates[posE,])

minuslogl<-function(b0=NULL,b1=NULL,b2=NULL,b3=NULL,b4=NULL,b5=NULL,b6=NULL,
b7=NULL,b8=NULL,b9=NULL, b10=NULL,b11=NULL, b12=NULL,b13=NULL,b14=NULL,b15=NULL,
b16=NULL,b17=NULL,b18=NULL,b19=NULL,
b20=NULL,b21=NULL,b22=NULL,b23=NULL,b24=NULL,b25=NULL,b26=NULL,
b27=NULL,b28=NULL,b29=NULL,b30=NULL,b31=NULL,b32=NULL,b33=NULL,
b34=NULL,b35=NULL,b36=NULL,b37=NULL,b38=NULL,b39=NULL,b40=NULL,
b41=NULL,b42=NULL,b43=NULL,b44=NULL,b45=NULL,b46=NULL,b47=NULL,
b48=NULL,b49=NULL,b50=NULL)
{
beta<-c(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11, b12,b13,b14,b15,b16,b17,b18,b19, 
b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b31, b33, b34, b35,
b36,b37,b38, b39, b40, b41, b42, b43, b44, b45, b46, b47, b48, b49, b50)

mllikpois <- -sum(covariatest%*%beta) + sum(exp(Bcovariates%*%beta ))
return(as.double(mllikpois))
}


n <- names(fixed)
fullcoef <- start
if (any(!n %in% names(fullcoef))) 
        stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
fullcoef[n] <- fixed
start[n] <- NULL
start <- sapply(start, eval.parent)
nm <- names(start)
oo <- match(nm, names(formals(minuslogl)) )
if (any(is.na(oo))) 
        stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
start <- start[order(oo)]
nm <- names(start)

f <- function(p) {
        l <- as.list(p)
        names(l) <- nm
        l[n] <- fixed
        do.call("minuslogl", l)
}
if (minfun=='nlminb')   oout <- nlminb(start=start, objective=f,...)
if (minfun=='optim')    oout <- optim(start, f, hessian = FALSE,...)
coef <- oout$par
names(coef)<-names(start)[order(oo)]
fullcoef[nm] <- coef
convergence<-oout$convergence
message<- oout$message
npar <- length(coef)
if (minfun=='optim') 
{
	min <- oout$value
        oout<-oout
	ooutb<- list()
}
if (minfun=='nlminb') 
{
	min <- oout$objective
        ooutb<-oout
	oout<- list()
}

fullcoef<-unlist(fullcoef)
if (is.null(n)) n<-'No  fixed parameters'
attr(fullcoef,'TypeCoeff')<-paste('Fixed:', n)


lambdafit<-as.numeric(exp(covariates%*%fullcoef))
 covariatesEst<-as.matrix(covariates[,order(oo)]) #only covariates with estimated parameters are kept in this matrix

if(modCI==TRUE) vcov<-VARbeta.fun(covariates=covariatesEst, lambdafit=lambdafit)
	else 	{vcov<-matrix(numeric(), 0L, 0L)
	 attr(vcov,'CalMethod')<-'Not required'}


if (modSim==FALSE)
{
cat(fill = TRUE)
cat('Number of observations  not used in the estimation process: ', (nobs-control),fill=TRUE)
cat('Total number of time observations: ',nobs,fill=TRUE)
cat("Number of events: ", length(posE), fill = TRUE)
cat(fill = TRUE)
cat('Convergence code: ',convergence, fill=TRUE)
if (convergence==0) cat('Convergence attained', fill=TRUE)
if (convergence!=0) 
{
cat(message, fill=TRUE)
stop('Convergence  not attained')
}
cat('Loglikelihood: ',-round(min,3), fill=TRUE)
cat(fill = TRUE)
cat('Estimated coefficients: ', fill=TRUE)
print(round(coef,3))
cat('Full coefficients: ', fill=TRUE)
print(round(fullcoef,3))
cat(fill = TRUE)
}

if (modCI==TRUE)
{

	if (CIty=='Transf')
	{
	CIlambda<-CItran.fun(VARbeta=vcov, lambdafit=lambdafit, covariates=covariatesEst,
	 clevel=clevel)
	}

	if (CIty=='Delta')
	{
	CIlambda<-CIdelta.fun(VARbeta=vcov, lambdafit=lambdafit, covariates=covariatesEst,
	 clevel=clevel)
	}
	UIlambda<-CIlambda$UIlambda
	LIlambda<-CIlambda$LIlambda
}
else
{
	UIlambda<-numeric(0)
	LIlambda<-numeric(0)
}

if (dplot==TRUE)
{
if ((modCI==TRUE)&(is.null(lambdaylim))) lambdaylim<-c(min(LIlambda, na.rm=TRUE), max(UIlambda, na.rm=TRUE))
plot(tim,lambdafit, ty='n',ylab='intensity', xlab='time', ylim=lambdaylim, xlim=lambdaxlim)
if (modCI==TRUE)
{
lines(tim,UIlambda, col='blue', lty=3)
lines(tim,LIlambda, col='red', lty=3)
legend(xlegend, legend=c('Fitted intensity', paste('Upper CI',CIty,sep=' '), paste('Lower CI',CIty,sep=' ')), col=c('black', 'blue', 'red'), lty=c(1,3,3), cex=0.8)
}
lines(tim,lambdafit, lwd=2)
mtext(paste(tit), outer = TRUE, line = -2,cex=1)
}

if (tind==TRUE)	covariates<-checkdim(covariates,j=1)

 new('mlePP',call = call, coef = coef, fullcoef = fullcoef, 
        vcov = vcov, min = min, details = oout, detailsb=ooutb, minuslogl = minuslogl, 
        nobs = if (missing(nobs)) 
            NA_integer_
        else nobs, method=minfun, npar=npar,inddat=inddat,
	lambdafit=lambdafit,  LIlambda=LIlambda, UIlambda=UIlambda,convergence=convergence,
	posE=posE, covariates=covariates, fixed=fixed, tit=tit,tind=tind, t=tim)

}





################################################################################
################################################################################
#         DENSITY ESTIMATION FOR NOISY OBSERVATIONS                            #
#            1) Known noise                                                    #
#            2) Unknown noise with sample of errors                            #
#            3) Unknown noise with sample of replicate noisy observations      #
#                                                                              #
# Julien Stirnemann                                                            #
# June 2012                                                                    #                  
################################################################################
################################################################################

################################################################################
#                        General functions                                     #
################################################################################

#define sinus cardinal
sinc  = function(x) {
y = matrix(1,dim(x)[1], dim(x)[2])
u = which(x!=0)
y[u] = sin(pi*x[u])/(pi*x[u])
return(y)
}

#define ifft (inverse FFT in Matlab)
ifft = function(x) fft(x, inverse=TRUE)/length(x)

#define function dlaplace
dlaplace <- function(x, mu=0, b=1) exp(-abs(x-mu)/b)/(2*b)

#define function rlaplace
rlaplace=function(n, mu=0, b=1) {
 u=runif(n, -0.5, 0.5)
 (mu-b*sign(u)*log(1-2*abs(u)))
 }


#calculation of the MISE
mise = function(density, obj){
if(class(obj)=="deamer"){
 supp <- obj$supp
 est <- obj$f
 } else { stop("argument 'obj' must be of class 'deamer'")} 

if(class(density)=="function"){
     dens <- density(supp)
     } else {
     stop("argument 'density' should be a function (see ?mise)")
   }

mise <- ((max(supp)-min(supp))/length(supp))*sum((dens-est)^2)
return(mise)
}


#positive part of the estimate
takepos=function(x) {x<-ifelse(x<0, 0, x); return(x)}

################################################################################
#                     Generic function methods                                 #
################################################################################

deamer <- function(y,errors, replicates, mu, sigma, noise.type, method,
                           grid.length=100, from, to, na.rm=FALSE) UseMethod("deamer")
plot.deamer <- function(x,type='l',main="",xlab='x',ylab='density',ylim,...) UseMethod("plot", x)
lines.deamer <- function(x,...) UseMethod("lines", x) 
print.deamer <- function(x,...) UseMethod("print", x)
predict.deamer <- function(object, newdata, na.rm=FALSE,...) UseMethod("predict", object)
          
          
deamer.default <- function(y, errors, replicates, mu, sigma, noise.type, method,
                           grid.length=100, from, to, na.rm=FALSE)
{
  x <- y
  if(!is.numeric(x))  stop("\targument 'y' must be numeric")
  x.name <- deparse(substitute(x))
  if(!is.vector(x)) warning("argument 'y' is coerced to a vector")
  x <- as.vector(x)
  if (any(is.na(x))) {
            if (na.rm) x <- x[!is.na(x)]
            else stop("argument 'y' contains missing values")
        }
  if (any(!is.finite(x))) {
   warning("non-finite values in 'y' have been removed")  
   x <- x[is.finite(x)]
  }
  n=length(x)
  
  if(missing(method)) stop("argument 'method' must be specified")  
  meth <- method
  if(missing(from)) from=min(x)
  if(missing(to))   to=max(x)
  supp <- seq(from, to, length=grid.length)
  
  if(meth%in%c("ke","se","ro")){
  if(meth=="ke"){
      if(missing(sigma)) stop("argument 'sigma' is missing with method='ke'") 
      if(!missing(errors)) warning("argument 'errors' is disregarded with method='ke'")
      if(!missing(replicates)) warning("argument 'replicates' is disregarded with method='ke'")
      M <- NULL; rep.name <- NULL 
      if(!is.numeric(sigma)) stop("sigma must be numeric")
      if(length(sigma)>1) stop("length(sigma) > 1")      
      if(sigma<=0) stop("sigma must be strictly positive")
      if(missing(mu)){
           warning("argument 'mu' is missing, default is mu=0") 
           mu <- 0
           }
      if(length(mu)>1) stop("length(mu) > 1")
      if(missing(noise.type)) {
         noise.type <- "Gaussian";
         warning("argument 'noise.type' is missing, default is 'Gaussian'")
         }      
      if(noise.type%in%c("Gaussian","gaussian","Laplace","laplace")==F) 
         stop("deamer does not support the specified error type")
      
      if(noise.type %in% c("Laplace","laplace")){
         invfestar <- make.invfestar.lap(n, mu, sigma)
         meth<-"kelap"
         } else {
            invfestar <- make.invfestar.gauss(n, mu, sigma)
            meth<-"kegauss"
            }
      Est <- est.ke(x, invfestar, supp)
   }
  
  if(meth=="se"){
       if(missing(errors)) stop("argument 'errors' is missing with method='se'") 
       rep.name=NULL;  
       if(any(!is.numeric(errors))) stop("argument 'errors' must be numeric")
       if(any(!is.finite(errors))) stop("argument 'errors' must be finite")
       if(!is.vector(errors)) warning("argument 'errors' is coerced to a vector")
       if(na.rm) {
          idx <- is.na(errors) 
          errors<-errors[-which(idx)]
          }
       M <- length(errors)
      if(M<2) stop("estimation is not possible with a sample of one error")
      if(!missing(mu)) 
            warning("argument 'mu' is disregarded with method='se'")
       if(!missing(sigma)) 
            warning("argument 'sigma' is disregarded with method='se'")
       if(!missing(noise.type)) 
            warning("argument 'noise.type' is disregarded with method='se'")
       if(!missing(replicates)) 
            warning("argument 'replicates' is disregarded with method='se'")
        sigma=NULL
        mu=NULL
       noise <- noise.char(eps=errors, n=n)
       Est <- est.se(x, noise, supp)
       }
             
  if(meth=="ro"){
     if(missing(replicates)) stop("argument 'replicates' is missing with method='ro'") 
     rep.name <- deparse(substitute(replicates))
     repl<-replicates
     repl <- as.matrix(repl)
     if(ncol(repl)< 2) stop("argument 'replicates' must be a 2-column matrix.")
     if(ncol(repl)> 2) stop("deamer does not yet support more than 2 replicates")
     if(sum(apply(repl,2,is.numeric)) < ncol(repl))  stop("argument 'replicates' must be numeric")      
     if (any(is.na(repl))) {
               if (na.rm){
                   allnona <- function(l) all(!is.na(l)) 
                   repl <- repl[apply(repl,1,allnona),]}
               else stop("argument 'replicates' contains missing values")
               }
     if (any(!is.finite(repl))) {
               allfinite <- function(l) all(is.finite(l))
               repl <- repl[apply(repl,1,allfinite),]
        }
     M=nrow(repl)
     if(M<2) stop("estimation is not possible with only one replicate")
      if(!missing(mu)) 
            warning("argument 'mu' is disregarded with method='ro'")
      if(!missing(sigma)) 
            warning("argument 'sigma' is disregarded with method='ro'")
      if(!missing(errors)) 
            warning("argument 'errors' is disregarded with method='ro'")
      if(!missing(noise.type)) 
            warning("argument 'noise.type' is disregarded with method='ro'")

      sigma=NULL; mu=NULL
      noise <- noise.char.repeat(repl[,1], repl[,2], x)
      Est <- est.ro(c(x,repl[,1]), noise, supp)
   }
  } else {stop("argument 'method' is not recognized as a deamer method")}
   
  phi <- sqrt(Est$m/pi)*sinc(Est$m/pi*matrix(1,511,1)%*%matrix(supp,nrow=1,ncol=length(supp))-
                                matrix(-255:255,511,1)%*%matrix(1,1,length(supp)))
  f <- Est$ahat%*%phi
  fout <- takepos(f)

  return(structure(list(y = x, f = fout,
                        n = n, M=M, method = meth,
                        mu=mu, sigma=sigma, 
                        supp = supp, m=Est$m, ahat=Est$ahat),
                        class="deamer")) 
                        
}

deamerKE <- function(y, mu, sigma, noise.type, grid.length=100, from, to, na.rm=FALSE){
  if(missing(noise.type)){
   noise.type <- "Gaussian"
   warning("argument 'noise.type' is missing, default is 'Gaussian'")
   }
  if(missing(sigma)) stop("argument 'sigma' is missing")
  if(missing(mu)){
           warning("argument 'mu' is missing, default is mu=0") 
           mu <- 0
           }
  deamer.default(y=y, mu=mu, sigma=sigma, noise.type=noise.type, method="ke", grid.length=grid.length, from=from, to=to, na.rm=na.rm)
 }

deamerSE <- function(y, errors, grid.length=100, from, to, na.rm=FALSE){
  if(missing(errors)) stop("argument 'errors' is missing")
  deamer.default(y=y, errors=errors, method="se", grid.length=grid.length, from=from, to=to, na.rm=na.rm)
}

deamerRO <- function(y, replicates, grid.length=100, from, to, na.rm=FALSE){
  if(missing(replicates)) stop("argument 'replicates' is missing")
  deamer.default(y=y, replicates=replicates, method="ro", grid.length=grid.length, from=from, to=to, na.rm=na.rm)
}

################################################################################
### Generic function for deamer class
################################################################################

# print function

print.deamer <- function(x,...){
  if (x$method == "kegauss"){ 
  cat("\nDeconvolution estimation with error ~ N(",x$mu,",",x$sigma,")",
      "\n\nData: ", x$n, " obs.\n\n", sep = "" )}
  if (x$method == "kelap"){ 
  cat("\nDeconvolution estimation with error ~ Lap(",x$mu,",",x$sigma,")",
      "\n\nData: ", x$n, " obs.\n\n", sep = "" )}
  if (x$method == "se"){ 
  cat("\nDeconvolution estimation with sample of errors",
      "\n\nData: ", x$n, " obs.; Errors: ",
      x$M, " obs.\n\n", sep = "" )}
  if (x$method == "ro"){ 
  cat("\nDeconvolution estimation with replicate observations",
      "\n\nData: ", x$n, " obs.; Replicates: ",
       x$M, " obs.\n\n",  sep = "")}     
}



# plot functions

plot.deamer <- function(x, type='l',main="",xlab='x',ylab='density',ylim,...){ 
  if(missing(ylim)) ylim=c(0,max(x$f)+(max(x$f)-min(x$f))/10)
  
  plot(x$supp, x$f, type=type, main=main, xlab=xlab, ylab=ylab, ylim=ylim,...)
  }

 
lines.deamer <- function(x,...){
  lines(x$supp,x$f,...)
}
  
# predict function to compute pdf(X=x0)

predict.deamer <- function(object, newdata, na.rm=FALSE,...){
 new.name <- deparse(substitute(newdata))
 if(any(!is.numeric(newdata))) stop("argument'newdata' must be numeric")
 if(any(is.na(newdata))) if(na.rm){
       newdata <- newdata[-which(is.na(newdata))] 
       warning("argument newdata has missing values")
       } else {
        stop("argument 'newdata' has missing values")
        }
 if(any(!is.finite(newdata))) stop("argument 'newdata' must be finite")
 obj.name <- deparse(substitute(object))
 newdata <- as.vector(newdata) 
 ahat <- object$ahat
 m <- object$m
 meth <- object$meth
 u <- length(newdata)
 xmat<-matrix(newdata,1,u)
 phi <- sqrt(m/pi)*sinc(m/pi*matrix(1,511,1)%*%xmat-
                                matrix(-255:255,511,1)%*%matrix(1,1,u))
 y <- ahat%*%phi
 return(as.vector(y))
}

################################################################################
################################################################################
#                      ROUTINES                                                #
################################################################################
################################################################################

################################################################################
#                   1) Known density of the noise                              #
################################################################################

#Ref: COMTE et LACOUR, "Data driven density estimation in presence of unknown 
#                              convolution operator"
#Ref: COMTE,TAUPIN and ROZENHOLC


#1.1) Compute the inverse of the characteristic function of the noise
#     (f^*_\varepsilon)^(-1)

## 1.1.a) Gaussian error
make.invfestar.gauss=function(n, mu, sigma){
mmax=100
invfestar.=matrix(NA, mmax, 256)
delta.=matrix(NA, mmax,1)

for(m in 1:mmax){
 lm=m/4
 z=matrix(lm*(2*(0:255)/256-1), 1,256)
 invfestar.[m,]=exp(-1i*mu*z+((sigma*z)^2)/2)
 delta.[m,]=2*lm/256*sum(abs(invfestar.[m,])^2)
}
mn=max(which(delta.<=n))
results=list(invfestar=invfestar., mn=mn, sigma=sigma, noisetype='gauss')
}

## 1.1.b) Laplace error
make.invfestar.lap=function(n, mu, sigma){
mmax=100
invfestar..=matrix(NA, mmax, 256)
delta..=matrix(NA, mmax,1)

for(m in 1:mmax){
 lm=m/4
 z=matrix(lm*(2*(0:255)/256-1), 1,256)
 invfestar..[m,]=((1+sigma^2*z^2)/exp(1i*mu*z))
 delta..[m,]=2*lm/256*sum(abs(invfestar..[m,])^2)
}
mn=max(which(delta..<=n))
results=list(invfestar=invfestar.., mn=mn, sigma=sigma, noisetype='lap')
}

#1.2) Density estimation
est.ke=function(Y, invfestar, supp){
l=4
  sigma=invfestar$sigma
  TYPE=invfestar$noisetype
  
  invqetoile=invfestar$invfestar
  mn=invfestar$mn
  u=length(supp)
  x=matrix(supp,1,u)
  n=length(Y)
	YY = matrix(Y,1,n)
  achapeau = matrix(0,mn,511)
  integrale=matrix(0,mn,1)
  gamman = matrix(0,mn,1)
  pen = matrix(0,mn,1)
  intfgauss=function(x, sigma) exp((sigma*x)^2)
  intflap=function(x, sigma) 1/(1+(sigma*x)^2)
  pengauss=function(x, LM, sigma) 0.5*LM^3*x
  penlap=function(x, LM, sigma) 4*(LM + 2/3/sigma^(-2)*LM^3+1/5/sigma^(-4)*LM^5)
  if(TYPE=='gauss') {functy=intfgauss; penty=pengauss}
  if(TYPE=='lap') {functy=intflap; penty=penlap}
   
	for (m in 1:mn){
		lm = m/l
    z = matrix(lm*(2*(0:255)/256-1),1,256)
		psiY = matrix(apply(exp(1i*t(z)%*%YY),1,sum),1,256)/n

    S=(psiY)*invqetoile[m,]
    achapeau[m,256:1]=sqrt(lm/pi)*Re(ifft(S))*(-1)^(0:255)
    achapeau[m,256:511]=sqrt(lm/pi)*Re(ifft(Conj(S)))*(-1)^(0:255)

    gamman[m,]=-sum(achapeau[m,]^2)
    integrale[m,]=sum(functy(lm*(0:99)/100, sigma))/100
    pen[m,]=penty(integrale[m,], lm, sigma)/n
    }

    mc=which.min(gamman[1:mn]+pen[1:mn])
    lmc=mc/4
#		phidex=sqrt(lmc/pi)*sinc(lmc/pi*matrix(1,511,1)%*%x-
#                                matrix(-255:255,511,1)%*%matrix(1,1,u))
#		fc=matrix(achapeau[mc,], 1,511)%*%phidex
		return(list(#f = as.vector(fc), 
    m=lmc, #supp = supp, 
    ahat=matrix(achapeau[mc,], 1,511)))
}

################################################################################
#           2) Unknown density of the noise, sample of pure errors             #
################################################################################

#
noise.char=function(eps, n){
M = length(eps);
Eps = matrix(eps,1,M)
mmax = 100
invqetoilechapeau = matrix(0,mmax,256)
invqetoiletilde = matrix(0,mmax,256)
deltatilde = matrix(0,mmax,1)

l=4

z.=matrix((2*(0:255)/256-1),1,256)
lm=matrix((1:mmax)/l,ncol=1)
z=lm%*%z.

for (l in 1:mmax){
    	invqetoilechapeau = (matrix(M/apply(exp(matrix(1i*z[l,], ncol=1)%*%Eps),1,sum),1,256))
    	invqetoiletilde[l,] = (invqetoilechapeau)*(abs(invqetoilechapeau)<M^0.5)
    	deltatilde[l,] = 2*lm[l,]/256*(sum(abs(invqetoilechapeau)^2))
}

mnchapeau=max(which(deltatilde<=(n/2)))
results=list(invqetoiletilde=invqetoiletilde, deltatilde=deltatilde, 
              mnchapeau=mnchapeau, mmax=mmax)
}


est.se=function(Y, noise, supp){
l=4
  invqetoiletilde=noise$invqetoiletilde
  deltatilde=noise$deltatilde
  mnchapeau=noise$mnchapeau
  mmax=noise$mmax
  u=length(supp)
  x=matrix(supp,1,u)
  n=length(Y)
	YY = matrix(Y,1,n)
  achapeau2 = matrix(0,mmax,511)
  gamman2 = matrix(0,mmax,1)
  pen2 = matrix(0,mmax,1)
	
	for (m in 1:mnchapeau){
		lm = m/l
    z = matrix(lm*(2*(0:255)/256-1),1,256)
		psiY = matrix(apply(exp(1i*t(z)%*%YY),1,sum),1,256)/n

    S2=(psiY)*invqetoiletilde[m,]
    achapeau2[m,256:1]=sqrt(lm/pi)*Re(ifft(S2))*(-1)^(0:255)
    achapeau2[m,256:511]=sqrt(lm/pi)*Re(ifft(Conj(S2)))*(-1)^(0:255)

    gamman2[m,]=-sum(achapeau2[m,]^2)
    pen2[m,]=2*(log(deltatilde[m,])/log(m+1))^2*deltatilde[m,]/n;
    }

		mc2=which.min(gamman2[2:mnchapeau]+pen2[2:mnchapeau]);
		lmc2=mc2/l
		phidex2=sqrt(lmc2/pi)*sinc(lmc2/pi*matrix(1,511,1)%*%x-
                                matrix(-255:255,511,1)%*%matrix(1,1,u))
		fc2=matrix(achapeau2[mc2,], 1,511)%*%phidex2
		return(list(f = as.vector(fc2), 
      m = lmc2, supp = supp, 
      ahat=matrix(achapeau2[mc2,], 1,511)))
  }



################################################################################
#                   3) Replicate noisy observations                            #
################################################################################

#DENSITY DECONVOLUTION WITH UNKNOWN NOISE AND AN AUXILIARY SET WITH 2 REPLICATES.
#PROJECTION WITH SINUS CARDINAL
#ADAPTATIVE METHOD

#Ref: COMTE et LACOUR, "Data driven density estimation in presence of unknown 
#                              convolution operator"
#Ref: COMTE, SAMSON and STIRNEMANN, HAL
#Ref: STIRNEMANN, COMTE and SAMSON, Statistics in Medicine 2012

#Y=X+e are independent noisy observations
#Y1 and Y2 are replicates of noisy observations


#Example:
# M=50 size of the replicate data set
# sig=0.5 SD for Gaussian noise
# Y1=rgamma(M , 12 , 0.8) + rnorm(M, 0, sigma); Y2=rgamma(M , 12 , 0.8) + 
#                                                            rnorm(M, 0, sigma)
# N=500 #size of non-replicate data-set
# Y=rgamma(N , 12 , 0.8) + rnorm(N, 0, sigma)
# supp=seq(0,40,length=100) #defines the support for estimation
# noise=noise.char.repeat(Y1,Y2) #estimation of empirical characteristic function of the noise
# f=estim.f(Y, noise, supp) #density estimation

# x11(); curve(dgamma(x , 12 , 0.8), from=0, to=40, lwd=2);
#        points(f~supp, type='l', lty=2)


#estimation of the inverse of the characteristic function of the noise (1/f^*_e)
#from the replicate dataset
noise.char.repeat=function(Y1, Y2, Y){
M = length(Y1)
n = length(Y);
eps = Y1-Y2
Eps = matrix(eps,1,M)
mmax = 100
invqetoilechapeau = matrix(0,mmax,256)
invqetoiletilde = matrix(0,mmax,256)
deltatilde = matrix(0,mmax,1)

l=4

z.=matrix((2*(0:255)/256-1),1,256)
lm=matrix((1:mmax)/l,ncol=1)
z=lm%*%z.

for (l in 1:mmax){
    	invqetoilechapeau = sqrt( abs( 
                            matrix(M/apply(cos(matrix(z[l,], ncol=1)%*%Eps),1,sum),1,256)
                            ))
    	invqetoiletilde[l,] = (invqetoilechapeau)*(abs((invqetoilechapeau))<M^0.25)
    	deltatilde[l,] = 2*lm[l,]/256*sum(abs(invqetoilechapeau)^2)
}

mnchapeau=max(which(deltatilde<=2*(n+M)))
results=list(invqetoiletilde=invqetoiletilde, deltatilde=deltatilde, 
              mnchapeau=mnchapeau, mmax=mmax)
}


# estimation of the density of X (named f) from Y+Y1

est.ro=function(Y, noise, supp){
l=4
  invqetoiletilde=noise$invqetoiletilde
  deltatilde=noise$deltatilde
  mnchapeau=noise$mnchapeau
  mmax=noise$mmax
  u=length(supp)
  x=matrix(supp,1,u)
  n=length(Y)
	YY = matrix(Y,1,n)
  achapeau2 = matrix(0,mmax,511)
  gamman2 = matrix(0,mmax,1)
  pen2 = matrix(0,mmax,1)
	
	for (m in 1:mnchapeau){
		lm = m/l
    z = matrix(lm*(2*(0:255)/256-1),1,256)
		psiY = matrix(apply(exp(1i*t(z)%*%YY),1,sum),1,256)/n

    S2=(psiY)*invqetoiletilde[m,]
    achapeau2[m,256:1]=sqrt(lm/pi)*Re(ifft(S2))*(-1)^(0:255)
    achapeau2[m,256:511]=sqrt(lm/pi)*Re(ifft(Conj(S2)))*(-1)^(0:255)

    gamman2[m,]=-sum(achapeau2[m,]^2)
    pen2[m,]=2*(log(deltatilde[m,])/log(m+1))^2*deltatilde[m,]/n;
    }

		mc2=which.min(gamman2[2:mnchapeau]+pen2[2:mnchapeau]);
		lmc2=mc2/l

		return(list( 
      m = lmc2,
      ahat=matrix(achapeau2[mc2,], 1,511)))
  }


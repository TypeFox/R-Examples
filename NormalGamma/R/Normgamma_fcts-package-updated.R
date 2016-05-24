require(optimx)	


dnormgam = function(par,x=NULL,N0=65536,plot=TRUE,log=FALSE, tail.cor=TRUE,cor=1e-15,mu=par[1],sigma=par[2],k=par[3],theta=par[4]) {
	
	# compute the density of X = S + B
	# where S ~ Gamma(k,theta) and B ~ N(mu,sigma^2)
	#
	# theta : scale of the gamma
	# k : shape of the gamma
	
	# The density is computed on [0,Tmax]
	
	
	if(sum(par!=c(mu,sigma,k,theta)))stop("Inconsistent definition of the parameters")
	if(sigma<0) stop("negative argument 'sigma'")
	if(k<0) stop("negative argument 'k'")
	if(theta<0) stop("negative argument 'theta'")


	xout <- x
	
	if (!is.null(xout)) {
		Trange = diff(range(xout))
		Tmin = min(xout)-0.1*Trange
		Tmax = max(xout)+0.1*Trange
	} else {
		Tmin = mu-5*sigma
		Tmax = mu+5*sigma+qgamma(0.99999,shape=k,scale=theta)
	}

	A = pi*N0/(Tmax-Tmin)
	
	mu = mu-Tmin
		
	L = -A +2*A*(0:(N0-1))/N0
		# L for lambda lives in frequency domain
	
	# characteristic function of the gamma S at L
	a = 1-1i*L*theta
	Cs = a^(-k)
	
	# characteristic function of the gaussian B at L
	b = -sigma^2*L^2/2+1i*mu*L
	Cb = exp(b)
	
	# characteristic function of X = S + B at L
	Y = Cs*Cb
	
	# Density of X = S + B at location T (time domain)
	T = pi*(0:(N0-1))/A
	
	fT = Re(A/N0/pi*exp(1i*A*T)*fft(Y))
	fT[fT<0] = 0

	# queue correction with gamma approximation
	if (tail.cor) {ind = which((T>mu+k*theta)&(abs(fT)<cor))
	}else{ind = NULL}
	
	if (length(ind)>0){
		ind = ind[1]:length(fT)
		badfT = fT[ind]
		fT[ind] = dgamma(T[ind]-mu,shape=k,scale=theta)
		}
	else badfT=NULL
	
	if (log){ # for log=TRUE only
		fT = log(fT)
		if (length(ind)>0) 
			fT[ind] = dgamma(T[ind]-mu,shape=k,scale=theta,log=TRUE)
	}
	
	T = T+Tmin
	jnd = which(T<Tmax)
	T = T[jnd]
	fT = fT[jnd]
		
	if (!is.null(xout)) {
		fT = approx(x=T,y=fT,xout=xout,rule=2)$y
		T = xout	
	}
	
	if (plot) plot(T,fT,type='l')

	list(dout=fT,xout=T)
	
	}
	
##############################################################

normgam.signal <- function(x,par,tail.cor=TRUE, cor=1e-15, gshift=FALSE,mu=par[1],sigma=par[2],k=par[3],theta=par[4]){
	
	# Estimation of Shat in the Normal-Gamma deconvolution
	# where X = S + N is observed with S~Gamma and N~Normal
	
	xout <- x
	Sh=xout # initialization
		
	tri=sort(xout,index.return=TRUE)
	xout2=tri$x
	Iout=tri$ix
	
	D1=dnormgam(par=c(par[1],par[2],par[3]+1,par[4]),x=xout2,tail.cor=tail.cor,cor=cor,plot=FALSE)
	D2=dnormgam(par=par,x=xout2,tail.cor=tail.cor,cor=cor,plot=FALSE)
	R = D1$dout/D2$dout

	# Due to numerical instability in close to 0, we have to adapt the left part of the curve
	# for very small relative values of xout 
	I0 = which.min(R); 
	if (I0>1)
		R[1:(I0-1)] = approx(c(0,I0),c(0,R[I0]),xout=1:(I0-1),rule=2)$y
	Shat=approx(D1$xout,k*theta*R,xout=xout2,rule=2)$y
	
	if ((gshift)&&(k<1)) { ### Not used in the paper
		print('Use gamma shift correction');
		# try to find the mode on [0,K*theta/10]
		jnd = which(Shat<k*theta/10)
		B = hist(Shat[jnd],breaks=5000,plot=FALSE)
		shift = B$mids[which.max(B$counts)]
		print(paste('shift =',shift))
		# shift from the mode (it is expected in 0)
		Shat <- Shat-shift
		Shat[Shat<0]=0
		}
	
	
	Sh[Iout] = Shat
	return(Sh)
	
	}
	
#############################################################



#############################################################

normgam.fit = function(X, N =NULL, par.init=NULL, lower = NULL, upper = NULL, control=NULL, verbose=FALSE){
	
Regular <- X
Negative <- N

## Is R version >= 3.0.0 ?
Rv3 <- (R.Version()$major>2)

mkdfftNormalGammaOptim = function(Regular,Negative) {
	
	# Construct the likelihood fct for the (Regular,negative) probes of a given array
	# using the transformed set of parameters
	
	n = length(Regular)+length(Negative)
	# Likelihood expressed in regular parameters
	function(part) {
		partinv = c(part[1],part[2],(part[3]/part[4])^2,part[4]^2/part[3])
		
		return(-(sum(dnormgam(partinv,x=Regular,log=TRUE,plot=FALSE)$dout)
		+sum(dnorm(Negative, mean=part[1],sd=part[2],log=TRUE)))/n )
	}
}

mkdfftNormalGammaOptim2 = function(Regular) {
	
	# Construct the likelihood fct for the Regular probes of a given array
	# using the transformed set of parameters
	
	n = length(Regular)
   # Likelihood expressed in regular parameters
	function(part) {
		partinv = c(part[1],part[2],(part[3]/part[4])^2,part[4]^2/part[3])
		return(-(sum(dnormgam(partinv,
			x=Regular,log=TRUE,plot=FALSE)$dout))/n)		}
}


# Define initial values of parameters

estimation_param_normexp_RMA <- function(pm , n.pts=2^14){
	
	set.seed(1234)
	max.density <- function(x, n.pts) {
		aux <- density(x, kernel = "epanechnikov", n = n.pts,
		na.rm = TRUE)
		aux$x[order(-aux$y)[1]]
		} # estimator of the mode of the density from sample x
	pmbg <- max.density(pm, n.pts)
	bg.data <- pm[pm < pmbg]
	pmbg <- max.density(bg.data, n.pts)
	bg.data <- pm[pm < pmbg]
	bg.data <- bg.data - pmbg
	bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2)
	sig.data <- pm[pm > pmbg]
	sig.data <- sig.data - pmbg

		alpha <- max.density(sig.data, n.pts)
	mubg <- pmbg	
	

	list(mu=mubg,sigma=bgsd,alpha=alpha)
	}


if(length(N)==0){
	dfftNormalGammaOptim = mkdfftNormalGammaOptim2(Regular)
	
	prma =estimation_param_normexp_RMA(Regular)	
    mu0=prma$mu
    sigma0 = prma$sigma
    theta0= ((sd(Regular))^2-sigma0^2)/(max(mean(Regular),median(Regular))-mu0)
    k0=(max(mean(Regular),median(Regular))-mu0)/theta0
    
    if(!((k0>0)&(theta0>0))){
     k0=1
     theta0= max(mean(Regular) - mu0,mu0/10) }


	}else{

dfftNormalGammaOptim = mkdfftNormalGammaOptim(Regular,Negative)
mu0 = median(Negative)
sigma0 = IQR(Negative)/1.349
theta0 = ((sd(Regular))^2-sigma0^2)/(mean(Regular)-mu0)
k0 = (mean(Regular)-mu0)/theta0 

if(!((k0>0)&(theta0>0))){
     k0=1
     theta0= max(mean(Regular) - mu0,mu0/10) }

          }
 
 
# transform parameter to be on data scale
ptrans <- function(U) c(U[1],U[2],U[3]*U[4],sqrt(U[3])*U[4])

#p3 = k0*theta0
#p4 = sqrt(k0)*theta0

if (verbose & (length(par.init)==0)) {
	print('Initial parameters mu0 sigma0 k0 theta0')
	print(paste(c(mu0,sigma0,k0,theta0),sep=','))
	}
	
if (verbose & (length(par.init)!=0)) {
	print('Initial parameters mu0 sigma0 k0 theta0')
	print(par.init)
	}	
	
	
	
# MLE optimization (MLEt with transformed parameters)

lower0 <- lower
upper0 <- upper
control0 <- control	
par.init0 <- par.init 

part0 = ptrans(c(mu0,sigma0,k0,theta0))

if(length(control)==0) control = list(maxit=1000,parscale=part0/10,dowarn=FALSE)



# Define initial new parameters


if( (length(par.init)==4)&(length(lower)==4)&(sum(lower<par.init) != 4)) stop('STOP : Initial values are out of bounds')
if( (length(par.init)==4)&(length(upper)==4)&(sum(upper>par.init) != 4)) stop('STOP : Initial values are out of bounds')

if(length(par.init)!=4){par.initt=part0
	}else{ par.initt = ptrans(par.init)}
	

		
if(length(N)==0){
	
		
   if(length(lower)==0) {lowert = c(-Inf,0,0,0) } else { lowert = ptrans(lower)}
   if(length(upper)==0) {uppert = c(max(Regular),sd(Regular)*3,part0[3]*3,part0[4]*3) } else  { uppert = ptrans(upper)}

   if(sum((lowert<par.initt)*(uppert>par.initt)) == 4 ){
          MLEt = optimx(par.initt, dfftNormalGammaOptim, method = "L-BFGS-B", lower = lowert, upper = uppert, control = control)
         if(Rv3){ coef.mlet <- coef(MLEt) } else { coef.mlet = unlist(MLEt$par)}
	}else{ coef.mlet=NA }
	
	if (sum(is.na(coef.mlet))){
	 if(length(lower0)+length(upper0) !=0){ print("Maximum out of bounds , Run maximization without box constraints")    }
	MLEt = optimx(par.initt,lower=c(-Inf,0,0,0), dfftNormalGammaOptim, method = "L-BFGS-B", control = control)
      if(Rv3){ coef.mlet <- coef(MLEt) } else { coef.mlet = unlist(MLEt$par)} 	
    }
	
	
	#if ((sum(is.na(coef.mlet))!=0)&(length(control0)>0)){
   # print('Maximization with default optimx control parameters ')
    
    
    	if (sum(is.na(coef.mlet))!=0){
  if((length(control0)>0)) print('Maximization with default optimx control parameters ')
MLEt = optimx(par.initt,lower=c(-Inf,0,0,0), dfftNormalGammaOptim, method = "L-BFGS-B" )
	if(Rv3){ coef.mlet <- coef(MLEt) } else { coef.mlet = unlist(MLEt$par)} 
	}
	
	if ((sum(is.na(coef.mlet))!=0)&(length(par.init)==4) ){
	  print('Convergence fails with initial values ->  Maximization with default initial values ')
	MLEt = optimx(part0,lower=c(-Inf,0,0,0), dfftNormalGammaOptim, method = "L-BFGS-B" )}


	

}else{
   	
    if(length(lower)==0) {lowert = c(min(Negative),sd(Negative)/3,part0[3]/3,part0[4]/3) }else{ lowert= ptrans(lower)}
    if(length(upper)==0) { uppert = c(max(Negative),sd(Negative)*3,part0[3]*3,part0[4]*3) } else { uppert = ptrans(upper)}
  
    if(sum((lowert<par.initt)*(uppert>par.initt)) == 4 ){	MLEt = optimx(par.initt, dfftNormalGammaOptim, method = "L-BFGS-B", lower = lowert,  upper = uppert, control = control)
	if(Rv3){ coef.mlet <- coef(MLEt) } else { coef.mlet = unlist(MLEt$par)} 
	
	     }else{coef.mlet=NA}
	
     
# if the optimization produces NA, restarts with non bounded optimization
if (sum(is.na(coef.mlet))){
	 if(length(lower0)+length(upper0) !=0){ print("Maximum out of bounds , Run maximization without box constraints")    }
	MLEt = optimx(par.initt,lower=c(-Inf,0,0,0), dfftNormalGammaOptim, method = "L-BFGS-B", control = control)
	
	if(Rv3){ coef.mlet <- coef(MLEt) } else { coef.mlet = unlist(MLEt$par)} 
    }

if ((sum(is.na(coef.mlet))!=0)&(length(control0)>0)){
    print('Maximization with default optimx control parameters ')
	MLEt = optimx(par.initt,lower=c(-Inf,0,0,0), dfftNormalGammaOptim, method = "L-BFGS-B" )
	if(Rv3){ coef.mlet <- coef(MLEt) } else { coef.mlet = unlist(MLEt$par)} 
	}

if ((sum(is.na(coef.mlet))!=0)&(length(par.init)==4) ){
	  print('Convergence fails with user initial values ->  Maximization with default initial values ')
	MLEt = optimx(part0,lower=c(-Inf,0,0,0), dfftNormalGammaOptim, method = "L-BFGS-B" )}



}           

# extract parameter estimates, likelihood and convergence code.
if(Rv3){ par <- coef(MLEt)
	     lik = MLEt$value
	     conv = MLEt$convcode
} else { par = unlist(MLEt$par,use.names=FALSE)
	 	 lik = as.numeric(MLEt$fvalues)
	 	 conv = as.numeric(MLEt$conv)
	 	 } 



# reverse transformation of parameters      
par <-  c(par[1],par[2], (par[3]/par[4])^2,par[4]^2/par[3] )

L=list()
L$par=par
L$lik=lik
L$conv=conv

return(L)		
		
}




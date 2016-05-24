#updated apr8, 2015

#-----------------------------------------
# regular cable fit deviance surface code
#-----------------------------------------

basic.cable<-function(z, g){

	if(g > 0)
		((z + g)^2/(4 * g)) * (z >=  - g & z <= g) + z * (z > g)
	else z * (z > 0)
}

fullcable.t<-function(t, b0, b1, b2, tau, gamma)
{
	b0 + b1 * t + b2 * basic.cable(t - tau, gamma)
}

cable.design.mat<-function(n,t.vect,tau,gamm)
{
	cbind(int=rep(1,n),t=t.vect,cable.t=basic.cable(t.vect-tau,gamm))
}


cable.loop<-function(i, j, tau.vect, gamm.vect, y.vect, design, max.flag = TRUE)
{
        t.hat <- tau.vect[i]
        g.hat <- gamm.vect[j]
        design[,3] <- basic.cable(design[,2] - t.hat, g.hat)
        fit <- try(lm(formula = y.vect ~ t + cable.t, data = design))

	#if(length(fit) == 12)  
   if(!inherits(fit,"try-error"))
       return( - dim(design)[1] * log(summary(fit)$sigma^2))
	else
		return(NA)
}

loglik.given.g<-function(i.vect, j, tau.vect, gamm.vect, y.vect,
	design, max.flag = TRUE){

		sapply(i.vect, cable.loop, j = j, tau.vect=tau.vect, gamm.vect=gamm.vect,
			y.vect=y.vect, design = design, max.flag = max.flag)
}

loglik.surf<-function(i.vect, j.vect, tau.vect, gamm.vect, y.vect,
	design, max.flag = TRUE){

	sapply(j.vect, loglik.given.g, i.vect = i.vect, tau.vect=tau.vect,
		gamm.vect=gamm.vect, y.vect=y.vect, design = design, max.flag = max.flag)
}

#---------------------------------------
# AR(p) cable fit deviance surface code
#---------------------------------------

is.stationary<-function(phi.vect){

	###############
	# stand-alone #
	########################################################
	# This function checks if the provided AR coefficients #
	# correspond to a stationary AR time series.           #
	########################################################

	#len<-length(try(arima.sim(n=2,list(ar=phi.vect)),silent=TRUE))

	#if(len==2)
   #   return(TRUE)
	#else
	#	return(FALSE)

    trial<-try(arima.sim(n=2,list(ar=phi.vect)),silent=TRUE)

    return(!inherits(trial,"try-error"))

}

ar.p.resid.core<-function(b0,b1,b2,tau,gamm,phi.vect,y.vect,t.vect,n,pp){

	if(pp==1)
		phi.vect<-as.matrix(phi.vect)

	err.vect<-y.vect-fullcable.t(t.vect,b0,b1,b2,tau,gamm)

	v.mat<-matrix(nrow=(n-pp),ncol=(pp+1))

	for(j in pp:0)
		v.mat[,(pp-j+1)]<-err.vect[(j+1):(n-pp+j)]
		
	innov<-as.vector(v.mat[,1]-v.mat[,-1]%*%phi.vect)

	return(list(resid=err.vect,innov=innov))

}

sse.ar.p.core<-function(b0,b1,b2,tau,gamm,phi.vect,y.vect,t.vect,n,pp){

	return(sum(
		ar.p.resid.core(b0,b1,b2,tau,gamm,phi.vect,y.vect,t.vect,n,pp)$innov^2
	))

}

sse.ar.p<-function(theta,y.vect,t.vect=NA,n=NA,pp=NA,stick=FALSE){

	# `optim()' doesn't like the variable `p'!

	b0<-theta[1]
	b1<-theta[2]
	b2<-theta[3]
	tau<-theta[4]

	if(!stick){

		gamm<-theta[5]
		phi.vect<-theta[-(1:5)]
	}

	else{

		gamm<-0
		phi.vect<-theta[-c(1:4)]
	}	

	if(is.na(n)){

		n<-length(y.vect)
		t.vect<-c(0:(n-1))

	}

	if(is.na(pp))
		pp<-length(phi.vect)

	return(sse.ar.p.core(b0,b1,b2,tau,gamm,phi.vect,y.vect,t.vect,n,pp))
}

sse.ar.p.given.beta.phi<-function(theta,beta.vect,phi.vect,y.vect,t.vect,n,pp){

	tau<-theta[1]
	gamm<-theta[2]
	b0<-beta.vect[1]
	b1<-beta.vect[2]
	b2<-beta.vect[3]

	return(sse.ar.p.core(b0,b1,b2,tau,gamm,phi.vect,y.vect,t.vect,n,pp))
}

sse.ar.p.given.t.g<-function(theta,tau,gamm,y.vect,t.vect,n,pp){

	b0<-theta[1]
	b1<-theta[2]
	b2<-theta[3]
	phi.vect<-theta[-(1:3)]

	return(sse.ar.p.core(b0,b1,b2,tau,gamm,phi.vect,y.vect,t.vect,n,pp))
}

sse.ar.p.given.g<-function(theta,gamm,y.vect,t.vect,n,pp){

	b0<-theta[1]
	b1<-theta[2]
	b2<-theta[3]
	tau<-theta[4]
	phi.vect<-theta[-(1:4)]

	return(sse.ar.p.core(b0,b1,b2,tau,gamm,phi.vect,y.vect,t.vect,n,pp))
}

cable.loop.ar.p<-function(tau,gamm,y.vect,t.vect,n,pp,fit=FALSE,verbose=FALSE)
{
	xreg<-cable.design.mat(n,t.vect,tau,gamm)[,-1]

	if(verbose)
		message(paste("tau=",tau,", gamma=",gamm))

	fitted<-arima(y.vect,order=c(pp,0,0),xreg=xreg,
		method="CSS",optim.control=list(maxit=5000))

	phi.hat<-fitted$coef[1:pp]
	beta.hat<-fitted$coef[-c(1:pp)]

	if(!fit){

		sse<-sse.ar.p.given.t.g(c(beta.hat,phi.hat),tau,gamm,y.vect,t.vect,n,pp)

		return(-(n-pp)*log(sse/(n-pp)))
	}

	else{

		return(fitted)

	}
}

loglik.given.g.ar.p<-function(tau.vect,gamm,y.vect,t.vect,n,pp){

	sapply(tau.vect,cable.loop.ar.p,gamm=gamm,
		y.vect=y.vect,t.vect=t.vect,n=n,pp=pp)
}

loglik.surf.ar.p<-function(
	tau.vect,gamma.vect,y.vect,t.vect=NA,n=NA,pp=NA){

	if(is.na(n)){

		n<-length(y.vect)
		t.vect<-c(0:(n-1))

	}

	sapply(gamma.vect,loglik.given.g.ar.p,
		tau.vect=tau.vect,y.vect=y.vect,t.vect=t.vect,n=n,pp=pp)

}

cable.dev<-function(tau.vect,gamma.vect,y.vect,t.vect=NULL,p=0){

	###############
	# stand-alone #
	####################################################################
	# 'tau.vect' and 'gamma.vect' specify a tau-gamma grid;            #
	# 'p' is the AR order for bent-cable residuals.                    #
	# This function computes the bent-cable profile deviance surface   #
	# matrix over the specified tau-gamma grid. It is intended to be   #
	# used for plotting a contour or perspective plot of the deviance  #
	# surface, in which case 'tau.vect' and 'gamma.vect' should have   #
	# length >=2 so that returned object is a matrix. If               #
	# such a plot is not required, then 'tau.vect' and 'gamma.vect'    #
	# can be scalar.                                                   #
	#                                                                  #
	# WARNING:                                                         #
	# When p>0, time series data are assumed, and 't.vect' should be   #
	# equidistant with unit increments.                                #
	####################################################################

	n<-length(y.vect)

	if(is.null(t.vect))
		t.vect<-c(0:(n-1))

	if(p==0){

		design<-data.frame(cable.design.mat(n,t.vect,0,0))

		loglik.mat<-loglik.surf(c(1:length(tau.vect)),c(1:length(gamma.vect)),
			tau.vect, gamma.vect, y.vect, design,TRUE)
	}

	else{

		message("Please be patient...")

		loglik.mat<-loglik.surf.ar.p(tau.vect,gamma.vect,y.vect,t.vect,n,p)
	}

	return(loglik.mat-max(na.omit(loglik.mat)))

}

#---------------------------
# plots code
#---------------------------

cable.ar.p.plot<-function(ar.p.fit,xlab="time",ylab="",
	main=NULL,ctp.ci=NULL){

	###############
	# stand-alone #
	#####################################################
	# This function plots the time series data for a    #
	# 'cable.ar.p.iter()' object provided by the user   #
	# as 'ar.p.fit'. This object must be obtained for   #
	# AR(p) data with p>0. Superimposed on the data is  #
	# the bent-cable fit from the object. The transition#
	# tau +/- gamma and tau are marked in red.          #
	# An optional 'cable.change.conf()' object as       #
	# 'ctp.ci' is also superimposed in blue.            #
	#####################################################


	y.vect<-ar.p.fit$y
	p<-ar.p.fit$p
	t.vect<-ar.p.fit$t

	main<-paste(main,": AR",p,"fit")

	plot(t.vect,y.vect,type="l",xlab=xlab,ylab=ylab,main=main)

	if(!ar.p.fit$stick)
		cable.lines(t.vect,ar.p.fit$est[1:5],col="red",lty=1)
	else
		cable.lines(t.vect,c(ar.p.fit$est[1:4],0),col="red",lty=1)

	if(!is.null(ctp.ci))
		abline(v=c(ctp.ci$change.hat,ctp.ci$interval),lty=2,col="blue")

}

cable.ar.p.diag<-function(ar.p.fit,resid.type="p",xlab="time",ylab="",
	main=NULL,main.all=NULL,ctp.ci=NULL){

	###############
	# stand-alone #
	#####################################################
	# This function produces (P)ACF diagnostic plots    #
	# for a 'cable.ar.p.iter()' object provided by the  #
	# user as 'ar.p.fit'. It plots the time series data #
	# from the object, and superimposes on it the bent- #
	# cable fit from the object. The transition tau +/- #
	# gamma and tau are marked in red.                  #
	# An optional 'cable.change.conf()' object as       #
	# 'ctp.ci' is also superimposed in blue. Plot title #
	# 'main' is for the time series plot, and           #
	# 'main.all' is for the entire set of plots.        #
	#####################################################

	par(mfcol=c(4,2))

	y.vect<-ar.p.fit$y
	t.vect<-ar.p.fit$t
	stick<-ar.p.fit$stick
	p<-ar.p.fit$p

	if(p==0){
		stop("This function cannot handle p=0.")
		return()
	}

	if(!stick)
		gamm<-ar.p.fit$est[5]
	else
		gamm<-0

	bend<-ar.p.fit$est[4]+c(-1,0,1)*gamm

	plot(1:2,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
	text(1.5,1.5,paste(main.all,"\n",main),cex=3)

	cable.ar.p.plot(ar.p.fit,xlab=xlab,ylab=ylab,main=main,ctp.ci=ctp.ci)

	fit.resid<-cable.ar.p.resid(ar.p.fit)

	plot(t.vect,fit.resid$resid,xlab=xlab,ylab="",
		main="fitted residuals",type=resid.type)
	abline(h=0,col="red")
	abline(v=bend,col="red")
	if(!is.null(ctp.ci))
		abline(v=c(ctp.ci$change.hat,ctp.ci$interval),lty=2,col="blue")

	plot(t.vect,c(rep(NA,p),fit.resid$innov),xlab=xlab,ylab="",
		main="fitted innovations",type=resid.type)
	abline(h=0,col="red")
	abline(v=bend,col="red")
	if(!is.null(ctp.ci))
		abline(v=c(ctp.ci$change.hat,ctp.ci$interval),lty=2,col="blue")

	acf(fit.resid$resid,main="ACF for Fitted Residuals")
	pacf(fit.resid$resid,main="PACF for Fitted Residuals")

	acf(fit.resid$innov,main="ACF for Fitted Innovations")
	pacf(fit.resid$innov,main="PACF for Fitted Innovations")

	par(mfrow=c(1,1))

}

cable.lines<-function(x,theta,col="black",lwd=1,lty=2,fit.lty=1){

	###############
	# stand-alone #
	####################################################
	# This function is meant to superimpose a fitted   #
	# bent cable over the plot of response-covariate   #
	# data the exhibit a bent cable structure.         #
	#                                                  #
	# 'x' contains the design points for the existing  #
	# plot, or the range of the design points.         #
	# This function then adds to it a bent cable with  #
	# coefficients 'theta'.                            #
	# Superimposed on the cable is the transition      #
	# marked by tau and tau +/- gamma. 'lwd' is for    #
	# the entire plot, 'lty' is for the transition,    #
	# and 'fit.lty' is for the cable.                  #
	####################################################

	changes<-theta[4]+c(-1,0,1)*theta[5]

	segments(min(x)-1,fullcable.t(min(x)-1,theta[1],theta[2],theta[3],
		theta[4],theta[5]),
		changes[1],fullcable.t(changes[1],theta[1],theta[2],theta[3],
		theta[4],theta[5]),col=col,lwd=lwd,lty=fit.lty)

	segments(changes[3],fullcable.t(changes[3],theta[1],theta[2],theta[3],
		theta[4],theta[5]),
		max(x)+1,fullcable.t(max(x)+1,theta[1],theta[2],theta[3],
		theta[4],theta[5]),col=col,lwd=lwd,lty=fit.lty)

	middle<-seq(changes[1],changes[3],range(changes)%*%c(-1,1)*.01)

	lines(middle,fullcable.t(middle,theta[1],theta[2],theta[3],
		theta[4],theta[5]),col=col,lwd=lwd,lty=fit.lty)

	abline(v=changes,col=col,lwd=lwd,lty=lty)

}


bentcable.dev.plot<-function(tau.vect,gamma.vect=NULL,y.vect,t.vect=NULL,stick=FALSE,p=0){

	##################
	# main interface #
	####################################################################
	# 'tau.vect' and 'gamma.vect' specify a tau-gamma grid.            #
	# 'p' is the AR order for bent-cable residuals.                    #
	# for 'stick=F': this function plots the contours of and returns   #
	#	the bent-cable profile deviance surface matrix over the         #
	#	specified tau-gamma grid - 'tau.vect' and 'gamma.vect'          #
	#	must have the same length                                       #
	# for 'stick=T': it plots and returns the broken-stick profile     #
	#	deviance function over the specified tau-domain while           #
	#	taking gamma=0 and thus ignores 'gamma.vect'                    #
	####################################################################

	n<-length(y.vect)

	if(is.null(t.vect))
		t.vect<-c(0:(n-1))

	if(!stick){

		dev<-cable.dev(tau.vect,gamma.vect,y.vect,t.vect,p)

		contour(tau.vect,gamma.vect,dev)
	}

	else{

		if(!is.null(gamma.vect))
			message("Ignoring 'gamma.vect'...")

		gamma.vect<-0

		dev<-cable.dev(tau.vect,gamma.vect,y.vect,t.vect,p)

		plot(tau.vect,dev,type="l",main="Broken-Stick Deviance",xlab="tau")

	}


	return(list(dev=dev,tau=tau.vect,gamma=gamma.vect))

}


cable.diagnose.resid<-function(fit){

	#plots PACF for residuals from "fit"
	#returns AR order based on AIC

	fit.resid<-resid(fit)

	pacf(fit.resid)

	return(c(p.aic=length(ar(fit.resid)$ar)))

}



#---------------------------
# AR(p) iterative fit code
#---------------------------

bentcable.ar<-function(y.vect,tgdev=NULL,p=0,stick=FALSE,t.vect=NULL,
	init.cable=NULL,init.phi=NULL,tol=1e-4,
	method0="css",method1="yw",ci.level=.95,main=NULL){

	##################
	# main interface #
	##################	
	
	# If 'init.cable' is missing then 'dev.mat'	
	# must be provided so as to produce an initial AR(p) fit.

	# 'method0' and 'method1' can be "css" (conditional sum-of-squares),
	# "mle", or "yw" (Yule-Walker). 'method0' is the default fitting method;
	# if 'method0' fails, then the fit is repeated using 'method1'.

	#####################################################################
	# This function fits a bent cable to response data 'y.vect' with    #
	# covariate values in 't.vect' (optional). For time-series data, a  #
	# stationary AR(p) autocorrelation structure is assumed. By         #
	# default, time points start at t=0, with increments of 1. Although #
	# the user can supply an alternative 't.vect', caution should be    #
	# taken in this case (see below).                                   #
	#                                                                   #
	# The fitted bent cable and estimated transition (tau and tau +/-   #
	# gamma) are plotted in red over the original response data. (P)ACF #
	# diagnostics are also plotted for fits with p>0. A confidence      #
	# interval for the critical time point (CTP) is also computed -     #
	# assuming the default time scale - and plotted in blue. By         #
	# definition, the CTP exists only when a change of sign in the      #
	# cable's slope takes place during transition. If the CTP (almost)  #
	# does not exist, then the CTP computation will fail. Failure can   #
	# also occur for a user-supplied time vector - even if successfully #
	# computed, the resulting CTP and confidence interval would be      #
	# meaningless as the computation assumes the default time scale.    #
	# The user is advised to rely on the default, and, if necessary,    #
	# transform the results to the preferred time scale after the model #
	# has been fitted.                                                  #
	#                                                                   #
	# Exception to this nuisance is in a non-time-series context where  #
	# the response data are independent, and the covariate is not time  #
	# and/or the design points are not equidistant. Then the "initial   #
	# AR(0) cable" is actually the max. likelihood fit for the original #
	# iid bent cable model (Chiu et al. 2006). No confidence intervals  #
	# are computed in this case.                                        #
	#                                                                   #
	# WARNING:                                                          #
	# The bent-cable likelihood surface often has multiple peaks. The   #
	# user is advised to examine the likelihood value to decide if the  #
	# returned fit corresponds to a local maximum. For a conditional    #
	# max. likelihood fit, the conditional sum-of-squares error (SSE)   #
	# can replace the likelihood for this purpose. If a 'css' fit is    #
	# successful, the conditional SSE is stored in                      #
	# '$cable$ar.p.fit$value' from the returned object. If 'css' fails  #
	# and the function refits using 'yw' or 'mle', then the conditional #
	# SSE is not directly retrievable. For model comparisons purposes,  #
	# the user can use the 'yw' or 'mle' estimate as the starting value #
	# for a subsequent 'css' fit, and the conditional SSE will appear   #
	# in the screen output as "initial value" while the 'css' algorithm #
	# iterates.                                                         #
	#####################################################################

	if(!is.null(tgdev) & (!is.null(init.cable) | !is.null(init.phi))){
		stop("'tgdev' cannot be used together with 'init.cable' or
			'init.phi'.")
		return()
	}

	if(p>0){

		if(is.null(init.phi))
			init.phi<-rep(c(.5,-.5),p)[1:p]

		else
			if(length(init.phi)!=p){

				message("****************************************")
				message("Your 'init.phi' and 'p' don't match:")
				message("ignoring 'init.phi', using given 'p'...")
				message("****************************************")

				init.phi<-rep(c(.5,-.5),p)[1:p]

			}

	}

	else

		if(length(init.phi)>0){

			message("****************************************")
			message("Your 'init.phi' and 'p' don't match:")
			message("ignoring 'p', using given 'init.phi'...")
			message("****************************************")

			p<-length(init.phi)

		}

	n<-length(y.vect)

	if(is.null(t.vect))
		t.vect<-c(0:(n-1))

	wrong.t<-FALSE

	if(sum(t.vect!=c(0:(n-1)))>0){

		message("*******************************")
		message("WARNING:")
		message("Cannot estimate CTP correctly")
		message("unless 't.vect' is c(0,1,2,...)")
		message("*******************************")

		wrong.t<-TRUE

	}

	if(!stick){ # fit bent cable:

		if(is.null(init.cable)){ # just getting initial values

			if(is.null(tgdev)){
				stop("Must supply 'init.cable' or 'tgdev'!")
				return()
			}

			message("*****************************************")
			message(paste("Finding initial values for AR(",p,") cable...",sep=""))
			message("*****************************************")

			if(p==0){

				message("*******************************")
				message("Finding best AR(0) cable fit...")
				message("*******************************")

				initfit<-try(cable.ar.0.fit(y.vect,t.vect,tgdev$tau,tgdev$gamma,tgdev$dev))

            #if(length(initfit$fit)!=6){
            if(inherits(initfit$fit,"try-error")){

					message("*****************************************************")
					message("AR(0) fit failed (deviance surface too irregular).")
					message("You may still try to use the printed values above")
					message("(if any) as initial values for subsequent AR(p) fits.")
					message("*****************************************************")
					return()

				}

				message("******************************************")
				message("If you have time-series data, then")
				message("choose your p according to AIC and PACF...")
				message("******************************************")

				print(cable.diagnose.resid(initfit$fit))

				message("************************************************")
				message("Please apply 'coef()[1:5]' to the '$fit' element")
				message("of the returned AR(0) fit as initial values for")
				message("subsequent AR(p) fits.")
				message("************************************************")

				return(initfit)

			}

			else{

				initfit<-try(cable.fit.known.change(y.vect,t.vect,n,tgdev$tau,tgdev$gamma,tgdev$dev,p))

				#if(length(initfit$fit)!=13){
            if(inherits(initfit$fit,"try-error")){

					message("*************************")
					message(paste("AR(",p,") fit failed.",sep=""))
					message("*************************")
					return()
				}

				message("**************************************")
				message("Please use the '$init' element of the")
				message("returned fit as initial values for")
				message(paste("subsequent AR(",p,") fits.",sep=""))
				message("**************************************")

				return(initfit)
				
			}
		}

		else{ # ready for real fit

			if(length(init.cable)>5){

				message("Using only 'init.cable[1:5]'...")
				init.cable<-init.cable[1:5]
			}

			init.vect<-c(init.cable,init.phi)

			fit<-cable.ar.p.iter(init.vect,y.vect,t.vect,n,tol,method0,method1)

			if(p==0){

				plot(t.vect,y.vect,xlab="",ylab="")
				title(paste("AR(0) bent cable:\n",main))
				cable.lines(t.vect,coef(fit$fit)[1:5],col="red")

			}

			if(!wrong.t)
				ctp.ci<-try(cable.change.conf(fit,ci.level))
			else
				ctp.ci<-NULL

			if(length(ctp.ci)!=3){

				message("******************************************")
				message("Failed to compute CTP confidence interval!")
				message("******************************************")

				if(p>0)
					cable.ar.p.diag(fit,main=main)

				return(list(cable=fit))
			}	

			else{

				if(p>0)
					cable.ar.p.diag(fit,main=main,ctp.ci=ctp.ci)
				else
					abline(v=c(ctp.ci$change,ctp.ci$int),lty=2,col="blue")

				return(list(cable=fit,ctp=ctp.ci))
			}
		}
	}

	# fit broken stick:


	if(is.null(init.cable)){ # just getting intial values

		if(is.null(tgdev)){
			stop("Must supply 'init.cable' or 'tgdev'!")
			return()
		}

		message("*****************************************")
		message(paste("Finding initial values for AR(",p,") stick...",sep=""))
		message("*****************************************")

		if(length(tgdev$gamma)>1 | tgdev$gamma[1]!=0){
			stop("Your 'tgdev' did not come from a stick fit.")
			return()
		}

		if(p==0){

			message("*******************************")
			message("Finding best AR(0) stick fit...")
			message("*******************************")

			initfit<-try(cable.ar.0.fit(y.vect,t.vect,tgdev$tau,tgdev$gamma,tgdev$dev,TRUE))
         
         #if(length(initfit$fit)!=6){
			if(inherits(initfit$fit,"try-error")){

				message("*****************************************************")
				message("AR(0) fit failed (deviance too irregular).")
				message("You may still try to use the printed values above")
				message("(if any) as initial values for subsequent AR(p) fits.")
				message("*****************************************************")
				return()

			}

			message("******************************************")
			message("If you have time-series data, then")
			message("choose your p according to AIC and PACF...")
			message("******************************************")

			message(cable.diagnose.resid(initfit$fit))

			message("************************************************")
			message("Please apply 'coef()[1:4]' to the '$fit' element")
			message("of the returned AR(0) fit as initial values for")
			message("subsequent AR(p) fits.")
			message("************************************************")

			return(initfit)

		}

		initfit<-try(cable.fit.known.change(y.vect,t.vect,n,tgdev$tau,tgdev$gamma,tgdev$dev,p))
      
      #if(length(initfit$fit)!=13){
		if(inherits(initfit$fit,"try-error")){

			message("*************************")
			message(paste("AR(",p,") fit failed.",sep=""))
			message("*************************")
			return()
		}

		message("**************************************")
		message("Please use the '$init' element of the")
		message("returned fit as initial values for")
		message(paste("subsequent AR(",p,") fits.",sep=""))
		message("**************************************")

		return(initfit)

	}

	# ready for real fit

	if(length(init.cable)>4){

		message("Using only 'init.cable[1:4]'...")
		init.cable<-init.cable[1:4]
	}

	if(p==0){

		fit<-stick.ar.0(init.cable,y.vect,t.vect,n)

		plot(t.vect,y.vect,xlab="",ylab="")
		cable.lines(t.vect,c(coef(fit$fit)[1:4],0),col="red")
		title(paste("AR(0) broken stick:\n",main))

	}

	else{	# fit AR broken stick:

		init.vect<-c(init.cable,init.phi)

		fit<-cable.ar.p.iter(init.vect,y.vect,t.vect,n,tol,method0,method1,stick=TRUE)
	}

	if(!wrong.t)
		ctp.ci<-try(cable.change.conf(fit,ci.level))
	else
		ctp.ci<-NULL

	if(length(ctp.ci)!=3){

		message("******************************************")
		message("Failed to compute CTP confidence interval!")
		message("******************************************")

		if(p>0)
			cable.ar.p.diag(fit,main=main)

		return(list(cable=fit))
	}	

	else{

		if(p>0)
			cable.ar.p.diag(fit,main=main,ctp.ci=ctp.ci)
		else
			abline(v=c(ctp.ci$change,ctp.ci$int),lty=2,col="blue")

		return(list(cable=fit,ctp=ctp.ci))
	}

}


cable.ar.0.fit<-function(y.vect,t.vect=NULL,tau.vect,gamma.vect,dev.mat,stick=FALSE){

	###############
	# stand-alone #
	#####################################################
	# This function fits a bent cable to response data  #
	# 'y.vect' over optional covariate value 't.vect',  #
	# assuming iid errors. Unlike 'bentcable.ar()',     #
	# this function does not require explicit initial   #
	# parameter values, but determines them via         #
	# 'tau.vect', 'gamma.vect', and 'dev.mat', which    #
	# should be obtained from 'bentcable.dev.plot()'.   #
   #####################################################


	n<-length(y.vect)

	if(is.null(t.vect))
		t.vect<-c(0:(n-1))

	if(!is.matrix(dev.mat))
		dev.mat<-matrix(dev.mat)

	fit0<-cable.fit.known.change(y.vect,t.vect,n,tau.vect,gamma.vect,dev.mat)

	if(!stick)
		return(cable.ar.p.iter(fit0$init,y.vect,t.vect,n,tol=NA))

	else
		return(stick.ar.0(fit0$init[1:4],y.vect,t.vect,n))

}


cable.fit.known.change<-function(y.vect,t.vect=NULL,n=NA,tau.vect,gamma.vect,dev.mat,p=0){

	###############
	# stand-alone #
	######################################################
	# This function fits a bent cable to response data   #
	# 'y.vect' over optional covariate values 't.vect',  #
	# conditioned on a known transition (i.e. known tau  #
	# and gamma). The user supplies 'tau.vect',          #
	# 'gamma.vect', and 'dev.mat' as a result of         #
	# 'bentcable.dev.plot()'. Using these, this function #
	# locates the peak of the deviance surface over the  #
	# specified tau-gamma grid, and computes the bent-   #
	# cable fit at this grid point. If multiple peaks    #
	# exist (such as along a ridge), then only that at   #
	# the smallest tau and smallest gamma is used.       #
	#                                                    #
	# For 'p=0', iid errors are assumed, and 't.vect'    #
	# can be non-equidistant. For p>0, AR(p) errors are  #
	# assumed, and 't.vect' must be equidistant with     #
	# unit increments, or the returned values will be    #
	# meaningless.                                       #
	######################################################


	peak<-pick.peak(dev.mat==0,nrow(dev.mat),ncol(dev.mat))

	t0<-tau.vect[peak$row]
	g0<-gamma.vect[peak$col]

	if(is.na(n))
		n<-length(y.vect)

	if(is.null(t.vect))
		t.vect<-c(0:(n-1))

	if(p==0){

		design<-data.frame(cable.design.mat(n,t.vect,t0,g0))

		fit<-lm(formula=y.vect~t+cable.t,data=design)

		return(list(fit=fit,init=c(coef(fit),t0,g0)))
	}

	else{

		fit<-cable.loop.ar.p(t0,g0,y.vect,t.vect,n,p,TRUE)

		init<-coef(fit)
		init<-c(init[-c(1:p)],t0,g0,init[1:p])

		names(init)[1:4]<-c("b0","b1","b2","tau")

		if(g0>0)
			names(init)[5]<-"gamm"
		else
			init<-init[-5]

		return(list(fit=fit,init=init))

	}

}

pick.peak<-function(dev.mat,nrow,ncol){

	pick.row<-dev.mat%*%rep(1,ncol)
	pick.row<-c(1:nrow)[pick.row>0][1]	# take smallest value in case of ridge

	pick.col<-rep(1,nrow)%*%dev.mat
	pick.col<-c(1:ncol)[pick.col>0][1]	# take smallest value in case of ridge

	return(list(row=pick.row,col=pick.col))
}


cable.ar.p.iter<-function(init,y.vect,t.vect=NULL,n=NA,
	tol,method0="css",method1="yw",stick=FALSE){

	###############
	# stand-alone #
	#########################################################################
	# for 'stick=F':                                                        #
	#	'init' must be in the order: b0, b1, b2, tau, gamma, phi.vect;       #
	#	if phi is present (i.e. AR(p) with p>0), then 't.vect' must be       #
	#	equidistant with unit increments                                     #
	#                                                                       #
	# for 'stick=T':                                                        #
	#	'init' must be in the order: b0, b1, b2, tau, phi.vect; phi          #
	#	must be present and 't.vect' must be equidistant with unit           #
	#       increments                                                      #
	#                                                                       #
	# 'method0' and 'method1' can be "css" (conditional sum-of-squares),    #
	#   "mle" or "yw" (Yule-Walker)                                         #
	#                                                                       #
	# 'tol' is the tolerance for determining convergence                    #
	#                                                                       #
	# This function fits a bent cable to the response data 'y.vect' with    #
	# optional covariate values in 't.vect', assuming AR errors. The AR     #
	# order is determined from the dimension of 'init'.                     #
	#########################################################################

	if(is.na(n))
		n<-length(y.vect)

	if(is.null(t.vect))
		t.vect<-c(0:(n-1))

	if(!stick)
		pp<-length(init[-c(1:5)])
	else
		pp<-length(init[-c(1:4)])

	if(pp==0){

		if(stick){
			stop("This function cannot handle AR(0) broken sticks.")
			return()
		}

		message("Trying 'nls()' Gauss-Newton algorithm...")

		fit0<-try(nls(formula=y.vect~fullcable.t(t.vect,b0,b1,b2,tau,gamm),trace=TRUE,
			start=list(b0=init[1],b1=init[2],b2=init[3],tau=init[4],gamm=init[5])))

      #if(length(fit0)!=6){
		if(inherits(fit0,"try-error")){

			message("Trying 'nls()' Golub-Pereyra algorithm...")

			fit0<-try(nls(formula=y.vect~fullcable.t(t.vect,b0,b1,b2,tau,gamm),trace=TRUE,
				algorithm="plinear",
				start=list(b0=init[1],b1=init[2],b2=init[3],tau=init[4],
				gamm=init[5])))

         #if(length(fit0)!=6){
			if(inherits(fit0,"try-error")){

				message("Trying 'nls()' Port algorithm...")

				fit0<-try(nls(formula=y.vect~fullcable.t(t.vect,b0,b1,b2,tau,gamm),trace=TRUE,
					algorithm="port",
					start=list(b0=init[1],b1=init[2],b2=init[3],tau=init[4],
					gamm=init[5])))

            #if(length(fit0)!=6){
            if(inherits(fit0,"try-error")){
					stop("Initial fit unsuccessful. Please try alternative initial values.")
					return()
				}
			}		
		}

		message("Converged!")

		return(list(fit=fit0,y=y.vect,t=t.vect,n=n,p=pp,stick=stick))

	}

	fit<-optim(init,sse.ar.p,y.vect=y.vect,t.vect=t.vect,n=n,pp=pp,stick=stick,
		control=list(trace=TRUE,abstol=tol,maxit=5000),method="BFGS")

	if(method0=="css"){

		if(!stick)
			phi.hat<-fit$par[-c(1:5)]
		else
			phi.hat<-fit$par[-c(1:4)]

		if(is.stationary(phi.hat)){

			message("CSS successful")
			return(list(estimate=fit$par,ar.p.fit=fit,y=y.vect,t=t.vect,n=n,p=pp,
				stick=stick,method=method0))
		}

		else{

			message("Non-stationary CSS fit")
			method0<-method1
		}

	}

	if(method0!="css"){

		message(paste("Applying method",method0,"..."))

		if(!stick)
			start<-fit$par
		else
			start<-c(fit$par[1:4],0,fit$par[-c(1:4)])

		b0<-start[1]
		b1<-start[2]
		b2<-start[3]
		tau<-start[4]
		gamm<-start[5]

		error<-999
		ct<-1

		while(error>tol){

			fit0.resid<-y.vect-fullcable.t(t.vect,b0,b1,b2,tau,gamm)

			Sig.ar.p.out<-cable.Sig.ar.p(fit0.resid,pp,n=n,method=method0)
			Sig<-Sig.ar.p.out$Sig
			phi.vect<-Sig.ar.p.out$phi

			xreg<-cable.design.mat(n,t.vect,tau,gamm)
			X.tran.Sig.inv<-t(xreg)%*%solve(Sig)

			beta.given.change<-solve(X.tran.Sig.inv%*%xreg)%*%
				X.tran.Sig.inv%*%y.vect

			b0<-beta.given.change[1]
			b1<-beta.given.change[2]
			b2<-beta.given.change[3]

			change.given.beta<-optim(c(tau,gamm),sse.ar.p.given.beta.phi,
				beta.vect=beta.given.change,phi.vect=phi.vect,
				y.vect=y.vect,t.vect=t.vect,n=n,pp=pp,
				control=list(trace=FALSE,abstol=tol,maxit=5000),method="BFGS")$par
				# conditional least squares

			tau<-change.given.beta[1]

			if(!stick)
				gamm<-change.given.beta[2]
			else
				gamm<-0

			if(gamm<0){

				message(paste("gamma hat =",gamm,"< 0: forcing to 2*tol"))
				gamm<-2*tol

			}


			est<-c(b0,b1,b2,tau,gamm,phi.vect)

			error<-sum(abs(start-est))
			message(paste(ct,": vector difference ",round(error,digits=ceiling(abs(log10(tol)))+2)))

			start<-est

			ct<-ct+1

		}

		names(est)[1:5]<-c("b0","b1","b2","tau","gamm")

		if(stick)
			est<-est[-5]

		return(list(estimate=est,ar.p.fit=Sig.ar.p.out$fit,y=y.vect,t=t.vect,n=n,p=pp,
			stick=stick,method=method0))
	}
}

cable.Sig.ar.p<-function(w.vect,p,n=NA,method="yw"){

	if(is.na(n))
		n<-length(w.vect)

	ar.p.fit<-ar(w.vect,aic=FALSE,order.max=p,demean=TRUE,method=method)
		# demean in case w.vect is from a fit that doesn't go through data

	#--------------------------------------------------------------
	# WARNING: ar() with method="mle" can give convergence problems
	# w/in iteration(s), although overall iterations converge
	#--------------------------------------------------------------

	phi.vect<-ar.p.fit$ar

	switch(method,

		yw={

			autocor.vect<-acf(w.vect,plot=FALSE,lag.max=p,demean=TRUE)$acf
		},

		mle={

			if(p==1)
				autocor.vect<-c(1,phi.vect[1])

			else{

				autocor.vect<-rep(0,p-1) # initialize from lags 1 to p-1
				Phi.mat<-matrix(nrow=(p-1),ncol=(p-1))
					# Phi.mat%*%autocor.vect = -phi.vect

				Phi.mat[1,1]<-phi.vect[2]-1
				if(p>2)
					Phi.mat[1,-1]<-phi.vect[-c(1:2)]

				phi.vect<-c(phi.vect,rep(0,p)) # 0's added for use in loop only

				for(h in 2:p-1){

					Phi.mat[h,h]<-phi.vect[2*h]-1
					Phi.mat[h,h-1]<-phi.vect[1]+phi.vect[2*h-1]

					if(h>2)
						Phi.mat[h,1:(h-2)]<-phi.vect[(h-1):2]+phi.vect[(h+1):(2*h-2)]

					if(2*h<p)
						Phi.mat[h,(h+1):(p-h)]<-phi.vect[(2*h+1):p]
				
				}

				Phi.mat<-replace(Phi.mat,is.na(Phi.mat),0)

				phi.vect<-phi.vect[1:p] # revert to original phi.vect

				autocor.vect<-solve(Phi.mat,-phi.vect[-p]) 
					# Phi.mat%*%autocor.vect = -phi.vect

				autocor.vect<-c(1,autocor.vect) # add back lag 0

				autocor.vect<-c(autocor.vect,phi.vect%*%autocor.vect[p:1]) 
					# add back lag p
			}

		}
	)

	for(i in (p+2):n)
		autocor.vect<-c(autocor.vect,phi.vect%*%autocor.vect[(i-1):(i-p)])

	Sig<-diag(n)

	for(i in 1:(n-1))
		Sig[i,((i+1):n)]<-autocor.vect[2:(n-i+1)]

	return(list(phi=phi.vect,Sig=(Sig+t(Sig)-diag(n)),fit=ar.p.fit))

}

cable.ar.p.resid<-function(ar.p.fit){

	###############
	# stand-alone #
	##################################################
	# This function computes fitted residuals and    #
	# innovations for a 'cable.ar.p.iter()' object   #
	# or the '$cable' element of a 'bentcable.ar()'  #
	# object, where p>0.                             #
	##################################################

	b0<-ar.p.fit$est[1]
	b1<-ar.p.fit$est[2]
	b2<-ar.p.fit$est[3]
	tau<-ar.p.fit$est[4]

	if(ar.p.fit$stick){

		phi.vect<-ar.p.fit$est[-c(1:4)]
		gamm<-0
	}

	else {

		gamm<-ar.p.fit$est[5]
		phi.vect<-ar.p.fit$est[-c(1:5)]
	}

	y.vect<-ar.p.fit$y
	t.vect<-ar.p.fit$t
	n<-ar.p.fit$n
	pp<-ar.p.fit$p

	return(ar.p.resid.core(b0,b1,b2,tau,gamm,phi.vect,y.vect,t.vect,n,pp))
}

stick.ar.0<-function(init.vect,y.vect,t.vect=NULL,n=NA){

	if(is.na(n))
		n<-length(y.vect)

	if(is.null(t.vect))
		t.vect<-c(0:(n-1))

	message("Trying 'nls()' Gauss-Newton algorithm...")

	fit<-try(nls(formula=y.vect~fullcable.t(t.vect,b0,b1,b2,tau,0),
		start=list(b0=init.vect[1],b1=init.vect[2],b2=init.vect[3],
		tau=init.vect[4]),trace=TRUE))
	
   #if(length(fit)!=6){
   if(inherits(fit,"try-error")){

		message("Trying 'nls()' Golub-Pereyra algorithm...")

		fit<-try(nls(formula=y.vect~fullcable.t(t.vect,b0,b1,b2,tau,0),
			algorithm="plinear",
			start=list(b0=init.vect[1],b1=init.vect[2],b2=init.vect[3],
			tau=init.vect[4]),trace=TRUE))

      #if(length(fit)!=6){
		if(inherits(fit,"try-error")){

			message("Trying 'nls()' Port algorithm...")

			fit<-try(nls(formula=y.vect~fullcable.t(t.vect,b0,b1,b2,tau,0),
				algorithm="port",
				start=list(b0=init.vect[1],b1=init.vect[2],b2=init.vect[3],
				tau=init.vect[4]),trace=TRUE))

         #if(length(fit)!=6){
			if(inherits(fit,"try-error")){
				stop("Initial fit unsuccessful. Please try alternative initial values.")
				return()
			}

		}
	}

	message("Converged!")

	return(list(fit=fit,y=y.vect,t=t.vect,n=n,p=0,stick=TRUE))
}


#----------------------------------------------------------------
# find confidence interval for point at which slope changes sign
#
# =================================================
# CANNOT HANDLE TIME VECTOR THAT'S NOT c(0,1,2,...)
# =================================================
#----------------------------------------------------------------

alpha<-function(theta.vect,summand,what){

	# computes alpha value in gradient formulas

	# `summand' is time point t
	# `what' is index of alpha

	what<-as.character(what)

	tau<-theta.vect[4]

	if(length(theta.vect)>4)
		gamm<-theta.vect[5]
	else
		gamm<-0

	return(
		switch(what,
			"1"=(summand>tau+gamm),
			"2"={
				if(gamm==0)
					0
				else
					(abs(summand-tau)<=gamm)*(summand-tau+gamm)/(2*gamm)
			},
			"3"={
				if(gamm==0)
					0
				else
					(abs(summand-tau)<=gamm)*(1-((summand-tau)/gamm)^2)/4
			},
			"4"=summand-tau
		)
	)
}

alph.vect<-function(theta.vect,p,summand,position){

	# produces alpha vector in gradient formulas

	# `summand' is time point t
	# `position' is variable wrt which derivative is taken
	# `position' can only be 3,4,5

	position<-as.character(position)

	b2<-theta.vect[3]

	if(length(theta.vect)>4)
		gamm<-theta.vect[5]
	else
		gamm<-0

	a.vect<-NULL

	switch(position,
		"3"={
			for(i in 0:p){
				a.vect<-c(a.vect,
					alpha(theta.vect,summand-i,1)*
					alpha(theta.vect,summand-i,4)+
					gamm*alpha(theta.vect,summand-i,2)^2)
			}	
		},

		"4"={
			for(i in 0:p){
				a.vect<-c(a.vect,
					alpha(theta.vect,summand-i,1)+
					alpha(theta.vect,summand-i,2))
			}	
		},

		"5"={
			for(i in 0:p){
				a.vect<-c(a.vect,
					alpha(theta.vect,summand-i,3))
			}
		}
	)

	return(a.vect)
}

grad<-function(theta.vect,phi.vect,summand,position){

	# takes gradient of score wrt theta

	# `summand' is time point t
	# `position' is variable wrt which derivative is taken

	position<-as.character(position)

	p<-length(phi.vect)

	phi.vect<-c(-1,phi.vect)
	b2<-theta.vect[3]

	return(
		switch(position,
			"1"=-sum(phi.vect),
			"2"=-sum(phi.vect*(summand-c(0:p))),
			"3"=-sum(phi.vect*alph.vect(theta.vect,p,summand,3)),
			"4"=b2*sum(phi.vect*alph.vect(theta.vect,p,summand,4)),
			"5"=-b2*sum(phi.vect*alph.vect(theta.vect,p,summand,5))
		)
	)

}

find.fisher<-function(n,theta.vect,phi.vect,sigma.sq){

	fisher.dim<-length(theta.vect)

	fisher<-matrix(0,ncol=fisher.dim,nrow=fisher.dim) # initialize
	working.fisher<-fisher

	for(summand in 0:(n-1)){
		for(j in 1:fisher.dim){
			for(k in j:fisher.dim){
				working.fisher[j,k]<-grad(theta.vect,phi.vect,summand,j)*
					grad(theta.vect,phi.vect,summand,k)
			}
		}
		fisher<-fisher+working.fisher
	}

	for(j in 2:fisher.dim)
		for(k in 1:(j-1))
			fisher[j,k]<-fisher[k,j]

	return(fisher/sigma.sq)
}



cable.change.conf<-function(ar.p.fit,level){

	###############
	# stand-alone #
	#####################################################
	# This function computes an approximate confidence  #
	# interval for the CTP of an AR(p) bent cable.      #
	#                                                   #
	# 'ar.p.fit' is a 'cable.ar.p.iter()' object, or    #
	# the '$cable' element of a 'bentcable.ar()'        #
	# object.                                           #
	#                                                   #
	# `level' is between 0 and 1.                       #
	#                                                   #
	# WARNING:                                          #
	# This function can only handle cases where time    #
	# steps are 0, 1, 2, ...                            #
	#####################################################

	stick<-ar.p.fit$stick

	if(!stick)
		k<-5
	else{
		k<-4
		gamm<-0
	}

	n<-ar.p.fit$n
	p<-ar.p.fit$p

	if(p>0){
	
		theta.vect<-ar.p.fit$est[1:k]
		phi.vect<-ar.p.fit$est[-c(1:k)]

		if(k==5)
			gamm<-theta.vect[5]

		if(ar.p.fit$method=="css")
			sigma.sq<-ar.p.fit$ar.p.fit$val/(n-p)
		else{
#			sigma.sq<-ar.p.fit$ar.p.fit$var
			sigma.sq<-sse.ar.p.core(theta.vect[1],theta.vect[2],theta.vect[3],
				theta.vect[4],gamm,phi.vect,ar.p.fit$y,ar.p.fit$t,n,p)/(n-p)
		}
	}

	else{

		theta.vect<-coef(ar.p.fit$fit)[1:k]
		phi.vect<-NULL

		fit.summ<-summary(ar.p.fit$fit)
		sigma.sq<-fit.summ$sig^2*fit.summ$df[2]/n

	}

	fisher.mat<-find.fisher(n,theta.vect,phi.vect,sigma.sq)

	tau<-theta.vect[4]
	
	if(!stick){

		b1<-theta.vect[2]
		b2<-theta.vect[3]
		gamm<-theta.vect[5]

		change.hat<-tau-gamm-2*b1*gamm/b2

		alpha.vect<-c(0,-2*gamm/b2,2*b1*gamm/b2^2,1,-(2*b1+b2)/b2)

		v<-alpha.vect%*%solve(fisher.mat)%*%alpha.vect
	}

	else{

		v<-solve(fisher.mat)[4,4]

		if(is.na(v)|v==0)
			return()

		change.hat=tau

	}	

	z<-qnorm(level+(1-level)/2)

	return(list(change.hat=as.numeric(change.hat),var=v,
		interval=change.hat+c(-1,1)*z*sqrt(v)))

}

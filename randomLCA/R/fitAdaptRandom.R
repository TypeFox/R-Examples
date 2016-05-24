`fitAdaptRandom` <- function(patterns,freq,nclass,calcSE,initoutcomep,initclassp,initlambdacoef,gh,constload,blocksize,probit,byclass,qniterations,penalty,verbose) {

# parameters
#   outcomes matrix of outcomes 0 or 1
#   freq vector of frequencies corresponding to each outcome combination
#   nclass number of classes
#   initoutcomep initial outcome probabilities
#   initclassp initial class probabilities
#   initlambdacoef initial lambda coefficient
#   constload are the loadings constant or different for each outcome
#   calcSE calculate standard errors ?
#   gh matrix of gauss-hermite coefficients first column positions, second columns weights
#   probit use probit transform rather than logitic to obtain outcome probabilities
#   verbose print information about algorithm    

	#browser()
	
	patterns <- as.matrix(patterns)
	mode(patterns) <- "double"

	nrepeats <- ifelse(constload,dim(patterns)[2],1)
	lambdasize <- ifelse(constload,1,min(dim(patterns)[2],blocksize))
	nlevel1 <- dim(patterns)[2]
	nlevel2 <- length(freq)

	if (verbose) print("fit.random.randomLCA")

		calclikelihood <- function(classx,outcomex,lambdacoef,momentdata,gh,patterns,calcfitted=FALSE,zprop=NULL) {
			#starttime <- proc.time()
			# turn classx into actual probabilities
			classp2 <- c(0,classx)       
			classp2 <- exp(classp2)/sum(exp(classp2))
			ill <- matrix(rep(NA,nclass*length(freq)),ncol=nclass)
			for (iclass in 1:nclass) {
				if (byclass) {
					#browser()
					ill[,iclass] <- .Call("bernoulliprobrandom",patterns,outcomex[iclass,],lambdacoef[iclass,],
						gh,momentdata,probit)*classp2[iclass]				} else {
					ill[,iclass] <- .Call("bernoulliprobrandom",patterns,outcomex[iclass,],lambdacoef,
						gh,momentdata,probit)*classp2[iclass]
				}
				##browser()
			}
# if zprop not supplied then we have the usual maximum likelihood
			if (is.null(zprop)) {
          ill2 <- rowSums(ill,na.rm=TRUE)
			    ll <- sum(log(ill2)*freq,na.rm=TRUE)
			} else {
# otherwise calculate the comple data maximum likelihood for the em algorithm
#        browser()
			  ill2 <- rowSums(log(ill)*zprop,na.rm=TRUE)
			  ll <- sum(ill2*freq,na.rm=TRUE)			  
			}
# penalise extreme outcome probabilities
			outcomep <- as.vector(1/(1+exp(abs(outcomex))))
#			pen <- dbeta(outcomep,1+penalty,1+penalty,log=TRUE)
			pen <-SciencesPo::ddirichlet(matrix(outcomep,nrow=1),rep(1+penalty/(nclass*2),length(outcomep)),log=TRUE)-SciencesPo::ddirichlet(matrix(outcomep,nrow=1),rep(1,length(outcomep)),log=TRUE)
			#			print(c(ll,sum(pen)))
			pen11 <- ll+sum(pen)
			if (is.nan(ll) || is.infinite(ll)) ll <- -1.0*.Machine$double.xmax
			if (calcfitted) {
# do this again in case we are using likelihood for em
			  ill2 <- rowSums(ill,na.rm=TRUE)
			  fitted <- ill2*sum(ifelse(apply(patterns,1,function(x) any(is.na(x))),0,freq))*
					ifelse(apply(patterns,1,function(x) any(is.na(x))),NA,1)
				classprob <- ill/ill2
				return(list(logLik=ll,penlogLik=pen11,fitted=fitted,classprob=classprob))
			} else return(list(logLik=ll,penlogLik=pen11))
		}  # end of calclikelihood

	calcrandom <- function(classx,outcomex,lambdacoef,momentdata) {

		classx <- c(0,classx)       
		classp <- exp(classx)/sum(exp(classx))
			
		onerandom <- function(x) {
	
				loglik <- function(beta,outcomes) {
			# calculate probabilities under each class
						for (i in 1:nclass) {
			# calculate the outcome probabilities for this class and current random
							if (byclass) {
								if (probit) outcomep <- pnorm(outcomex[i,]+rep(lambdacoef[i,],nrepeats)*beta)
								else outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef[i,],nrepeats)*beta))				
							} else {
								if (probit) outcomep <- pnorm(outcomex[i,]+rep(lambdacoef,nrepeats)*beta)
								else outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef,nrepeats)*beta))
							}
							oneprob <- exp(sum(outcomes*log(outcomep)+(1-outcomes)*log(1-outcomep),na.rm=TRUE))
			# multiply by class probabilities
							if (i==1) allprob <- oneprob*classp[i]
							else allprob <- allprob+oneprob*classp[i]
						}
					ll <- -(sum(log(allprob))+dnorm(beta,mean=0,sd=1,log=TRUE))
					if (is.nan(ll) || is.infinite(ll)) ll <- .Machine$double.xmax
				  return(ll)
				}
			  optim.fit <- nlm(loglik,x[length(x)],print.level=0,iterlim=1000,hessian=TRUE,outcomes=x[1:(length(x)-1)],gradtol=1.0e-7)
#  calculate se
			  return(c(beta=optim.fit$estimate[1],sebeta=sqrt(1/optim.fit$hessian)))
		}
		betas <- t(apply(cbind(patterns,momentdata[,1]),1,onerandom))
		return(betas)
	}
	
	adaptivefit <- function(classx,outcomex,lambdacoef,calcSE,momentdata,gh,patterns) {
	
	
		fitparams <- function(classx,outcomex,lambdacoef,
			momentdata,calcSE,gh,patterns,noiterations=qniterations,zprop=NULL) {
			calcllfornlm <- function(params,momentdata,gh,patterns,zprop) {
				oneiteration <- calclikelihood(if (nclass==1) NULL else params[1:(nclass-1)],
					matrix(params[nclass:(nclass+nlevel1*nclass-1)],nrow=nclass),
					if (byclass) matrix(params[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*lambdasize-1)],nrow=nclass)
					else params[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+lambdasize-1)],
					momentdata,gh,patterns,zprop=zprop)
				ll <- -oneiteration$penlogLik
				if (is.nan(ll) || is.infinite(ll)) ll <- .Machine$double.xmax
				ll
			}
			
			nlm1 <- nlm(calcllfornlm, c(classx, as.vector(outcomex), lambdacoef), iterlim = noiterations,
				print.level=ifelse(verbose,2,0),hessian=calcSE,stepmax=1,
				check.analyticals = FALSE,momentdata=momentdata,gh=gh,patterns=patterns,zprop=zprop)
			return(list(penlogLik=-nlm1$minimum,
				classx=(if (nclass==1) NULL else nlm1$estimate[1:(nclass-1)]),
				outcomex=matrix(nlm1$estimate[nclass:(nclass+nlevel1*nclass-1)],nrow=nclass),
				 lambdacoef=(if (byclass) matrix(nlm1$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*lambdasize-1)],nrow=nclass)
				else lambdacoef=nlm1$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+lambdasize-1)]),
				nlm=nlm1))	
		}
	
		#browser()
		
		oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns)
		currll <- oneiteration$penlogLik
	# shift the quadrature points for the first time
		nshift <- 0
		repeat {
		  momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
		  oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns)
		  # check if moving quadrature points has changed likelihood
		  if (verbose) cat("current ll",oneiteration$penlogLik,"\n")
		  if (abs((oneiteration$penlogLik-currll)/oneiteration$penlogLik) < 1.0e-6) {
		    currll <- oneiteration$penlogLik
		    break()
		  }
		  currll <- oneiteration$penlogLik
		  nshift <- nshift+1
		  if (nshift > 200) stop("too many shift iterations - increase quadrature points")
		}

		
		oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns,calcfitted=TRUE)
		currll <- oneiteration$penlogLik
    zprop <- oneiteration$classprobs

      
		adaptive <- TRUE
		prevll <- -Inf
		nadaptive <- 0
    while(adaptive) {
			# need to do an optimisation on the other parameters
			fitresults <- fitparams(classx,outcomex,lambdacoef,momentdata,FALSE,gh,patterns,zprop=zprop)
			currll <- fitresults$penlogLik
 #     print(cat("em likelihood",currll)
			outcomex <- fitresults$outcomex
			classx <- fitresults$classx
			lambdacoef <- fitresults$lambdacoef
			if (verbose) cat("current ll from optimisation",currll,"\n")
			optll <- currll
			#browser()
			# shift the quadrature points again
			nshift <- 0
			repeat {
				momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
				oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns)
				# check if moving quadrature points has changed likelihood
				if (verbose) cat("current ll",oneiteration$penlogLik,"\n")
				if (abs((oneiteration$penlogLik-currll)/oneiteration$penlogLik) < 1.0e-6) {
					currll <- oneiteration$penlogLik
					break()
				}
				currll <- oneiteration$penlogLik
				nshift <- nshift+1
       			if (nshift > 200) stop("too many shift iterations - increase quadrature points")
			}
        	adaptive <-	(abs((currll-prevll)/currll)>1.0e-7) || (abs((currll-optll)/currll)>1.0e-7)
        	if ((prevll-currll)/abs(currll) > 1.0e-3) stop("divergence - increase quadrature points")
        	nadaptive <- nadaptive+1
        	if (nadaptive > 500) stop("too many adaptive iterations - increase quadrature points")
        	prevll <- currll
    # get the proportions for the em algorithm
			oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns,calcfitted=TRUE)
      zprop <- oneiteration$classprobs
    }
		fitresults <- fitparams(classx,outcomex,lambdacoef,momentdata,calcSE,gh,patterns,noiterations=500)
		return(list(nlm=fitresults$nlm,momentdata=momentdata))
	} # end adaptivefit

	# momentdata is level2
	# mu2,sigma2,
	
	
    if (nclass==1) classx <- NULL
    else  {
        classx <- rep(NA,nclass-1)
    	initclassp <- ifelse(initclassp==0.0,1.0e-3,initclassp)       	
    	initclassp <- ifelse(initclassp==1.0,1-1.0e-3,initclassp)  
    	initclassp <- initclassp/sum(initclassp)
        for (i in 2:nclass) classx[i-1] <- log(initclassp[i]/initclassp[1])
    }

	initoutcomep <- ifelse(initoutcomep<1.0e-2,1.0e-2,initoutcomep)
	initoutcomep <- ifelse(initoutcomep>(1-1.0e-2),1-1.0e-2,initoutcomep)
    if (probit) outcomex <- qnorm(initoutcomep)
    else outcomex <- log(initoutcomep/(1-initoutcomep))


	momentdata <- matrix(rep(c(0,1),each=nlevel2),nrow=nlevel2)

	if (missing(initlambdacoef) || is.null(initlambdacoef)) lambdacoef <- rep(3.0,lambdasize)
	else lambdacoef <- initlambdacoef
	#browser()  	   
 	myfit <- adaptivefit(classx,outcomex,lambdacoef,calcSE,momentdata,gh,patterns)

	#browser()
	
	optim.fit <- myfit$nlm
	momentdata <- myfit$momentdata

	classx <- NULL
	if (nclass>1) classx <- optim.fit$estimate[1:(nclass-1)]
	outcomex <- matrix(optim.fit$estimate[nclass:(nclass+nclass*nlevel1-1)],ncol=dim(patterns)[2])
# transform using logistic to probabilities     
       if (byclass) lambdacoef <- matrix(optim.fit$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*lambdasize-1)],nrow=nclass)
       else lambdacoef <- optim.fit$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+lambdasize-1)]
	
	final <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns,calcfitted=TRUE)
			

# calculate the probabilities
# add extra column to classx
	classx <- c(0,classx)       
	 classp <- exp(classx)/apply(matrix(exp(classx),nrow=1),1,sum)

    if (probit) outcomep <- pnorm(outcomex)
    else outcomep <- exp(outcomex)/(1+exp(outcomex))
	
# extract the se
	if (!calcSE) separ <- rep(NA,length(optim.fit$estimate))
	else {
		s <- svd(optim.fit$hessian)
		separ <- sqrt(diag(s$v %*% diag(1/s$d) %*% t(s$u)))
		separ[!is.finite(separ)] <- NA
	}

# determine random effects

	ranef <- calcrandom(classx,outcomex,lambdacoef,momentdata)
	
    np <- length(optim.fit$estimate)
    nobs <- sum(freq)
    deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/final$fitted)))
    list(fit=optim.fit,nclass=nclass,classp=classp,outcomep=outcomep,lambdacoef=lambdacoef,se=separ,
    	np=np,nobs=nobs,logLik=final$logLik,penlogLik=final$penlogLik,observed=freq,fitted=final$fitted,
    	deviance=deviance,ranef=ranef,classprob=final$classprob)
}

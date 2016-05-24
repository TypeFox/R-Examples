
setClass("summary.bild", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric",call="language"))

setClass("bild", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric", residuals="numeric", 
		s.residuals="numeric",ind.probability="numeric", prob.matrix="matrix", Fitted="numeric", 
		Fitted.av="numeric", Time="numeric", model.matrix= "matrix", y.matrix="matrix",
		subset.data="data.frame", y.av="numeric", f.value="factor",call="language"))


setGeneric("getAIC",def=function(object) standardGeneric("getAIC"))


setGeneric("getLogLik",def=function(object) standardGeneric("getLogLik"))


bild<-function(formula = formula(data), data, time,id, subSET, aggregate=FALSE, start = NULL, trace = FALSE,
dependence="ind", method="BFGS", control=bildControl(), 
integrate=bildIntegrate())
{
# *****************DEFINITION OF INTERNAL FUNCTIONS ******************
# na.action for binary families

na.discrete.replace1 <- function(frame,  n.times, ti.repl)
	{
	vars <- names(frame)
	names(vars) <- vars
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	badlines<-NULL
	bad.ind<-NULL
	for(j in 1:length(vars)) 
	{k1<-1
	for (i in 1:n.cases)
	{k2<-cumti.repl[i]
	x <- frame[[j]][k1:k2]
	pos <- is.na(x)
	if(any(pos))
			if(j == 1){distance.between.na <- diff(seq(1, n.times)[!pos])
					if (any(distance.between.na > 2 ))
					{badlines<-c(badlines,c(k1:k2))
					bad.ind<-c(bad.ind,i)
					x[pos] <- -1}
				x[pos] <- -1}
			else stop("NA's on covariates not allowed")
	frame[[j]][k1:k2]<-x
	k1<-k2+1
	} 
	}
		if (length(bad.ind)>=1)
		cat("Warning Message: Condition on NA's not respected:\nresults might be inaccurate\n")
		return(list(data=frame, badlines=badlines,bad.ind=bad.ind))
	}




na.discrete.replace2 <- function(frame,  n.times, ti.repl)
	{
	vars <- names(frame)
	names(vars) <- vars
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	badlines<-NULL
	bad.ind<-NULL
	for(j in 1:length(vars)) 
	{k1<-1
	for (i in 1:n.cases)
	{k2<-cumti.repl[i]
	x <- frame[[j]][k1:k2]
	pos <- is.na(x)
	if(any(pos))
			if(j == 1){distance.between.na <- diff(seq(1, n.times)[!pos])
				if (any(distance.between.na > 2 ))
				{badlines<-c(badlines,c(k1:k2))
				bad.ind<-c(bad.ind,i)
				x[pos] <- -1}
				if (any(distance.between.na == 2 ) & length(distance.between.na)>1)
				{for (i2 in 1:(length(distance.between.na)-1))
				{if (distance.between.na[i2]==2 & distance.between.na[i2+1]==2)
				{badlines<-c(badlines,c(k1:k2))
				bad.ind<-c(bad.ind,i)
				x[pos] <- -1}}}
			x[pos] <- -1}
			else stop("NA's on covariates not allowed")
	frame[[j]][k1:k2]<-x
	k1<-k2+1
	} 
	}
		if (length(bad.ind)>=1)
		cat("Warning Message: Condition on NA's not respected:\nresults might be inaccurate\n")
		return(list(data=frame, badlines=badlines,bad.ind=bad.ind))
	}

#
# compute the fitted values for variuos families of distributions
#
	various.fitted <- function(Fitted)
	{
		#return(1/(1 + care.exp( - Fitted)))
		return(1/(1 +( - Fitted)))
	}
#
# compute numerical second derivatives
#
	num.info <- function(coefficients, FUN, X, data)
	{
		FUN <- get(FUN, inherits = TRUE)
		values <- FUN(coefficients, X, data)
		p <- length(values)
		Info <- matrix(0, p, p)
		h <- rep(0, p)
		delta <- cbind((abs(coefficients) + 1e-012) * 0.0001, rep(1e-012, p))
		delta <- apply(delta, 1, max)
		for(i in 1:p) {
			h[i] <- delta[i]
			new.values <- FUN(coefficients + h, X, data)
			Info[, i] <- (new.values - values)/delta[i]
			h[i] <- 0
		}
		Info
	}

# compute numerical second derivatives involving integrals

	num.infoI <- function(coefficients, FUN, X, data,integrate)
	{
		FUN <- get(FUN, inherits = TRUE)
		values <- FUN(coefficients, X, data,integrate)
		p <- length(values)
		Info <- matrix(0, p, p)
		h <- rep(0, p)
		delta <- cbind((abs(coefficients) + 1e-012) * 0.0001, rep(1e-012, p))
		delta <- apply(delta, 1, max)
		for(i in 1:p) {
			h[i] <- delta[i]
			new.values <- FUN(coefficients + h, X, data,integrate)
			Info[, i] <- (new.values - values)/delta[i]
			h[i] <- 0
		}
		Info
	}


########################################################################
# ------------------------------------------------------------------- 
	
logL.bin0.aux<- function(parameters, X, data, trace)
{
	loglik1<- function(param, X, y, trace)
	{#calculate logLi(beta/bi) for each individual
	npar <- as.integer(length(param)+1)
	beta <- as.double(param[1:(npar-1)])
	log.psi <- as.double(0)
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
	x <-matrix(as.double(X),nrow=n,ncol=(npar-1))
	theta <- work <- pij<-as.double(rep(0,n))
	logL <- L<-as.double(0)

		results <- .Fortran("mblik1",logL,pij,beta,log.psi, npar,x,y,theta,work,n,PACKAGE="bild")

	
	return(list(loglik=results[[1]],pij=results[[2]]))
	}
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	prob<-as.double(rep(0,length(y)))
	counts<-data[[3]]
	logL<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
	logL<-logL+counts[i]*z$loglik
	prob[k1:k2]<-z$pij
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")

return(list(nloglik=-logL, pij=prob))}
############

logL.bin0<- function(parameters, X, data,trace)
        logL.bin0.aux(parameters, X, data,trace)$nloglik
# ------------------------------------------------------------------- 
	
logL.bin1mf.aux<- function(parameters, X, data, trace)
{
	loglik1<- function(param, X, y, trace)
	{#calculate logLi(beta/bi) for each individual
	npar <- as.integer(length(param))
	beta <- as.double(param[1:(npar-1)])
	log.psi <- as.double(param[npar])
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
	x <-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta <- work <- pij<-as.double(rep(0,n))
	logL <- L<-as.double(0)

		results <- .Fortran("mblik1",logL,pij,beta,log.psi, npar,x,y,theta,work,n,PACKAGE="bild")

	return(list(loglik=results[[1]],pij=results[[2]]))
	}
		if(trace)	cat(paste(format(parameters[length(parameters)], digit=4), collapse=" "), "\t")
	
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	prob<-as.double(rep(0,length(y)))
	counts<-data[[3]]
	logL<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
	logL<-logL+counts[i]*z$loglik
	prob[k1:k2]<-z$pij
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")

return(list(nloglik=-logL, pij=prob))}
############

logL.bin1mf <- function(parameters, X, data,trace)
        logL.bin1mf.aux(parameters, X, data,trace)$nloglik

# ------------------------------------------------------------------- 

logL.bin2mf.aux <- function(parameters, X, data, trace)
{
	loglik2 <- function(param, X, y,trace)
	{
	npar <-as.integer(length(param))
	beta <- as.double(param[1:(npar-2)])
	lpsi <- as.double(param[(npar-1):npar])
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
	x <-matrix(as.double(X),nrow=n,ncol=npar-2)
	theta <- work <- pij<-as.double(rep(0,n))
	logL <- L<-as.double(0)

		results <- .Fortran("blik2m",logL,pij,beta,lpsi,npar,x,y,theta,work,n,PACKAGE="bild")


	return(list(loglik=results[[1]],pij=results[[2]]))}

		if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=3), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=3)), collapse=" "), "\t")

	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	prob<-as.double(rep(0,length(y)))
	counts<-data[[3]]
	logL<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- loglik2(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
	logL<-logL+counts[i]*z$loglik
	prob[k1:k2]<-z$pij
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")

return(list(nloglik=-logL, pij=prob))}

############

logL.bin2mf<- function(parameters, X, data,trace)
        logL.bin2mf.aux(parameters, X, data,trace)$nloglik

######################################################################## começa novo
# ------------------------------------------------------------------- 
	
logL.bin0Ifm<- function(parameters, X, data, integrate, trace)
{
	loglik1 <- function(param, X, y,integrate,trace)
	{
	npar <-as.integer(length(param))
	beta<- as.double(param[1:(npar-1)])
	bt<- as.double(param[1:(npar-1)])
	log.psi<-as.double(0)
	omega<-as.double(param[npar])
	y[is.na(y)]<-(-1)
	y<- as.integer(y)
	n <- as.integer(length(y)) 
	x<-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta<- work<-as.double(rep(0,n))
	logL <- as.double(0)

	li<-as.double(integrate$li)
	ls<-as.double(integrate$ls)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

 		results <- .Fortran("integ1",logL,bt,beta,log.psi,omega,npar,
		x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")

	return(results[[1]])}
	
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

	nparam<-length(parameters)
	omega1<-as.double(parameters[nparam])
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	mt<-as.double(rep(0,length(y)))
	logL1<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate,trace=trace)
	#logL1 gives the log-likelihood
	logL1<-logL1+counts[i]*z
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
return(nloglik=-logL1)}


######## to compute pij

logL.bin0I.aux<- function(parameters, X, data, trace)
{
	loglik1 <- function(param, X, y,trace)
	{
	
	npar <-as.integer(length(param))
	beta<- as.double(param[1:(npar-1)])
	bt<- as.double(param[1:(npar-1)])
	log.psi<-as.double(0)
	omega<-as.double(param[npar])
	y<- as.integer(y)
	n <- as.integer(length(y)) 
	x<-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta<- work<-pij<- prob<-as.double(rep(0,n))
	logL <- as.double(0)
	eta<-i.fit<- fit<-as.vector(npar-1)
	eta<-x%*%beta
	m<-glm(as.numeric(y)~offset(eta), family=binomial)
	bi<-coef(m)
	beta[1]<-beta[1]+bi
	y[is.na(y)]<-(-1)
	y<- as.integer(y)
	
	results <- .Fortran("mblik1",logL,prob,beta,log.psi, npar,x,y,theta,work,n,PACKAGE="bild")


	i.fit<-x%*%beta

	return(list(loglik=results[[1]],pij=results[[2]], fit=i.fit))}
	
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

	nparam<-length(parameters)
	omega1<-as.double(parameters[nparam])
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	mt<-fitted<-as.double(rep(0,length(y)))
	logL1<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],trace=trace)
	mt[k1:k2]<-z$pij
	fitted[k1:k2]<-z$fit
	k1<-k2+1
	}

return(list(pij=mt, fit=fitted))}


# ------------------------------------------------------------------- acaba novo


logL.bin1Ifm<- function(parameters, X, data, integrate, trace)
{
	loglik1 <- function(param, X, y,integrate,trace)
	{
	npar <-as.integer(length(param)-1)
	beta<- as.double(param[1:(npar-1)])
	bt<- as.double(param[1:(npar-1)])
	log.psi<-as.double(param[npar])
	omega<-as.double(param[npar+1])
	y[is.na(y)]<-(-1)
	y<- as.integer(y)
	n <- as.integer(length(y)) 
	x<-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta<- work<-as.double(rep(0,n))
	logL <- as.double(0)

	li<-as.double(integrate$li)
	ls<-as.double(integrate$ls)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

 		results <- .Fortran("integ1",logL,bt,beta,log.psi,omega,npar,
		x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")

	return(results[[1]])}
	
		if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=4), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

	nparam<-length(parameters)
	omega1<-as.double(parameters[nparam])
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	mt<-as.double(rep(0,length(y)))
	logL1<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate,trace=trace)
	#logL1 gives the log-likelihood
	logL1<-logL1+counts[i]*z
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
return(nloglik=-logL1)}


######## to compute pij

logL.bin1.aux<- function(parameters, X, data, trace)
{
	loglik1 <- function(param, X, y,trace)
	{
	
	npar <-as.integer(length(param)-1)
	beta<- as.double(param[1:(npar-1)])
	bt<- as.double(param[1:(npar-1)])
	log.psi<-as.double(param[npar])
	omega<-as.double(param[npar+1])
	y<- as.integer(y)
	n <- as.integer(length(y)) 
	x<-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta<- work<-pij<- prob<-as.double(rep(0,n))
	logL <- as.double(0)
	eta<-i.fit<- fit<-as.vector(npar-1)
	eta<-x%*%beta
	m<-glm(as.numeric(y)~offset(eta), family=binomial)
	bi<-coef(m)
	beta[1]<-beta[1]+bi
	y[is.na(y)]<-(-1)
	y<- as.integer(y)
	
	results <- .Fortran("mblik1",logL,prob,beta,log.psi, npar,x,y,theta,work,n,PACKAGE="bild")


	i.fit<-x%*%beta

	return(list(loglik=results[[1]],pij=results[[2]], fit=i.fit))}
	
		if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=4), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

	nparam<-length(parameters)
	omega1<-as.double(parameters[nparam])
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	mt<-fitted<-as.double(rep(0,length(y)))
	logL1<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],trace=trace)
	mt[k1:k2]<-z$pij
	fitted[k1:k2]<-z$fit
	k1<-k2+1
	}

return(list(pij=mt, fit=fitted))}

# ------------------------------------------------------------------- 

logL.bin2Ifm<- function(parameters, X, data,integrate,trace)
{
	loglik2 <- function(param, X, y,integrate,trace)
	{
	npar <-as.integer(length(param)-1)
	beta<- as.double(param[1:(npar-2)])
	bt<- as.double(param[1:(npar-2)])
	log.psi<-as.double(param[(npar-1):npar])
	omega<-as.double(param[npar+1])
	y[is.na(y)]<-(-1)
	y<- as.integer(y)
	n <- as.integer(length(y)) 
	x<-matrix(as.double(X),nrow=n,ncol=npar-2)
	theta<- work<- as.double(rep(0,n))
	logL <- as.double(0)

	li<-as.double(integrate$li)
	ls<-as.double(integrate$ls)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

 		results <- .Fortran("integ",logL,bt,beta,log.psi,omega,npar,
		x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")


	return(results[[1]])}

		if(trace)	cat(paste(format(parameters[length(parameters)-2], digit=4), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)-1], digit=4)), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

	nparam<-length(parameters)
	omega1<-as.double(parameters[nparam])
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	mt<-as.double(rep(0,length(y)))
	logL1<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik2(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate,trace=trace)
	#logL1 gives the log-likelihood
	logL1<-logL1+counts[i]*z
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")

return(nloglik=-logL1)}


########to compute pij

logL.bin2.aux<- function(parameters, X, data,trace)
{
	loglik2 <- function(param, X, y,trace)
	{
	npar <-as.integer(length(param)-1)
	beta<- as.double(param[1:(npar-2)])
	bt<- as.double(param[1:(npar-2)])
	lpsi<-as.double(param[(npar-1):npar])
	omega<-as.double(param[npar+1])
	y<- as.integer(y)
	n <- as.integer(length(y)) 
	x<-matrix(as.double(X),nrow=n,ncol=npar-2)
	theta<- work<- pij<-prob<-as.double(rep(0,n))
	logL <- as.double(0)
	eta<-i.fit<-fit<-as.vector(npar-2)
	eta<-x%*%beta
	m<-glm(as.numeric(y)~offset(eta), family=binomial)
	bi<-coef(m)
	beta[1]<-beta[1]+bi
	y[is.na(y)]<-(-1)
	y<- as.integer(y)

		results <- .Fortran("blik2m",logL,pij,beta,lpsi,npar,x,y,theta,work,n,PACKAGE="bild")

	i.fit<-x%*%beta

	return(list(loglik=results[[1]],pij=results[[2]],fit=i.fit ))}

		if(trace)	cat(paste(format(parameters[length(parameters)-2], digit=4), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)-1], digit=4)), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	prob<-fitted<-as.double(rep(0,length(y)))
	counts<-data[[3]]
	logL<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- loglik2(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
	prob[k1:k2]<-z$pij
	fitted[k1:k2]<-z$fit
	k1<-k2+1
	}

return(list(pij=prob, fit=fitted))}

########################################################################
# compute gradient loglik for binary response
#
gradlogL.bin0<- function(parameters, X,data, trace)
{
	gradient1<- function(param, X, y)
	{
	npar <- as.integer(length(param)+1)
	beta <- as.double(param[1:(npar-1)])
	log.psi <- as.double(0)
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
	theta <- work <- as.double(rep(0,n))
	gbeta <- dbeta <- dbeta1<- as.double(rep(0,npar-1))
	glpsi<- as.double(0)
	db <- matrix(as.double(0),nrow=3,ncol=npar-1)
	x <- matrix(as.double(X),nrow=n, ncol=npar-1)
	der<-as.double(rep(0,npar-1))

  		results <- .Fortran("mbgd1",gbeta,glpsi,beta,log.psi, 
             			npar,x,y,theta,work,der,db,dbeta,dbeta1,n,PACKAGE="bild")

		
	return(results[[1]])}

	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	gradlogL<-0
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- gradient1(param=parameters, X=X[k1:k2,],y=y[k1:k2])
	gradlogL<-gradlogL+counts[i]*z
	k1<-k2+1
	}
return(-gradlogL)}	

# ------------------------------------------------------------------- 
gradlogL.bin1mf<- function(parameters, X,data, trace)
{
	gradient1<- function(param, X, y)
	{
	npar <- as.integer(length(param))
	beta <- as.double(param[1:(npar-1)])
	log.psi <- as.double(param[npar])
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
	theta <- work <- as.double(rep(0,n))
	gbeta <- dbeta <- dbeta1<- as.double(rep(0,npar-1))
	glpsi<- as.double(0)
	db <- matrix(as.double(0),nrow=3,ncol=npar-1)
	x <- matrix(as.double(X),nrow=n, ncol=npar-1)
	der<-as.double(rep(0,npar-1))

  		results <- .Fortran("mbgd1",gbeta,glpsi,beta,log.psi, 
             			npar,x,y,theta,work,der,db,dbeta,dbeta1,n,PACKAGE="bild")

		
	gradL<-c(results[[1]],results[[2]])
	return(gradL)
	}
	nparam <- as.integer(length(parameters))
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	gradlogL<-0
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- gradient1(param=parameters, X=X[k1:k2,],y=y[k1:k2])
	gradlogL<-gradlogL+counts[i]*z
	k1<-k2+1
	}
return(-gradlogL)}	

# ------------------------------------------------------------------- 

gradlogL.bin2mf<- function(parameters, X,data, trace)
{
	gradient2 <-  function(param,X,y)
	{
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
  	npar <- as.integer(length(param))
	beta <- as.double(param[1:(npar-2)])
	lpsi <- as.double(param[(npar-1):npar])
  	theta <- work <- as.double(rep(0,n))
  	g.beta <- d.beta <- d.beta1<-d.beta2 <- der <- as.double(rep(0,npar-2))
  	g.lpsi1<-g.lpsi2 <- as.double(0)
  	db <- matrix(as.double(0),3,npar-2)
  	db1 <- matrix(as.double(0),4,npar-2)
  	db2 <- matrix(as.double(0),5,npar-2)
 	x <- matrix(as.double(X),nrow=n, ncol=npar-2)
  	
	result <- .Fortran("bgd2m",g.beta,g.lpsi1,g.lpsi2,beta,
	lpsi,npar,x,y,theta,work,d.beta,d.beta1,d.beta2,n,der,db,db1,db2,PACKAGE="bild")

	return(c(result[[1]],result[[2]],result[[3]]))
	}
	nparam <- as.integer(length(parameters))
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	gradlogL<-0
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- gradient2(param=parameters, X=X[k1:k2,],y=y[k1:k2])
	gradlogL<-gradlogL+counts[i]*z
	k1<-k2+1
	}
return(-gradlogL)}
# ------------------------------------------------------------------- 


######################################################################## começa novo
# ------------------------------------------------------------------- 

gradlogL.bin0Ifm <- function(parameters, X,data,integrate,trace)
{
	gradient1 <-  function(param,X,y,integrate)
	{
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
  	npar <- as.integer(length(param))
	beta <- as.double(param[1:(npar-1)])
	bt <- as.double(param[1:(npar-1)])
	lpsi <- as.double(0)
	omega<-as.double(param[npar])
  	theta <- work <- as.double(rep(0,n))
  	g.beta <- d.beta <- d.beta1<- der <- as.double(rep(0,npar-1))
  	g.lpsi1<- as.double(0)
  	d.beta <- matrix(as.double(0),3,npar-1)
	gvar<-as.double(0)
 	x <- matrix(as.double(X),nrow=n, ncol=npar-1)
  	db <- matrix(as.double(0),3,npar-1)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)


	result <- .Fortran("gint1",g.beta,g.lpsi1,gvar,bt,beta,lpsi,omega,npar,x,y,theta,work,n,
                  d.beta,d.beta1,der,db,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")


	return(c(result[[1]],result[[3]]))}

	loglik1 <- function(param, X, y,integrate)
	{
	npar <- as.integer(length(param))
	beta <- as.double(param[1:(npar-1)])
	bt <- as.double(param[1:(npar-1)])
	lpsi <- as.double(0)
	omega<-as.double(param[npar])
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
	x <-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta <- work <-as.double(rep(0,n))
	logL <- as.double(0)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

 		results <- .Fortran("integ1",logL,bt,beta,lpsi,omega,npar,
		x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")

return(results[[1]])}

	nparam <- as.integer(length(parameters))
	omega1<-parameters[nparam]
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	dbeta<-as.double(rep(0,nparam-1))
	dlog.psi1<-0
	dvar<-0
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
	grad<-gradient1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)

	for (j in 1:(nparam-1))
	{
	dbeta[j]<-dbeta[j]+counts[i]*(grad[j]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
	}

 #using the chain rule
	dvar<-dvar+counts[i]*(grad[nparam]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))*exp(omega1)
	k1<-k2+1
	}
gr<-c(dbeta, dvar)
 return(-gr)}


######################################################################## acaba novo
# ------------------------------------------------------------------- 


gradlogL.binI1fm <- function(parameters, X,data,integrate,trace)
{
	gradient1 <-  function(param,X,y,integrate)
	{
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
  	npar <- as.integer(length(param)-1)
	beta <- as.double(param[1:(npar-1)])
	bt <- as.double(param[1:(npar-1)])
	lpsi <- as.double(param[npar])
	omega<-as.double(param[npar+1])
  	theta <- work <- as.double(rep(0,n))
  	g.beta <- d.beta <- d.beta1<- der <- as.double(rep(0,npar-1))
  	g.lpsi1<-g.lpsi2 <- as.double(0)
  	d.beta <- matrix(as.double(0),3,npar-1)
	gvar<-as.double(0)
 	x <- matrix(as.double(X),nrow=n, ncol=npar-1)
  	db <- matrix(as.double(0),3,npar-1)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)


	result <- .Fortran("gint1",g.beta,g.lpsi1,gvar,bt,beta,lpsi,omega,npar,x,y,theta,work,n,
                  d.beta,d.beta1,der,db,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")


	return(c(result[[1]],result[[2]],result[[3]]))}

	loglik1 <- function(param, X, y,integrate)
	{
	npar <- as.integer(length(param)-1)
	beta <- as.double(param[1:(npar-1)])
	bt <- as.double(param[1:(npar-1)])
	lpsi <- as.double(param[npar])
	omega<-as.double(param[npar+1])
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
	x <-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta <- work <-as.double(rep(0,n))
	logL <- as.double(0)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

 		results <- .Fortran("integ1",logL,bt,beta,lpsi,omega,npar,
		x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")

return(results[[1]])}

	nparam <- as.integer(length(parameters)-1)
	omega1<-parameters[nparam+1]
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	dbeta<-as.double(rep(0,nparam-1))
	dlog.psi1<-0
	dvar<-0
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
	grad<-gradient1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)

	for (j in 1:(nparam-1))
	{
	dbeta[j]<-dbeta[j]+counts[i]*(grad[j]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
	}
	dlog.psi1<-dlog.psi1+counts[i]*(grad[nparam]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))		
 #using the chain rule
	dvar<-dvar+counts[i]*(grad[nparam+1]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))*exp(omega1)
	k1<-k2+1
	}
gr<-c(dbeta,dlog.psi1,dvar)
 return(-gr)}

# ------------------------------------------------------------------- 

gradlogL.binI2fm <- function(parameters, X,data,integrate,trace)
{
	gradient2 <-  function(param,X,y,integrate)
	{
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
  	npar <- as.integer(length(param)-1)
	beta <- as.double(param[1:(npar-2)])
	bt <- as.double(param[1:(npar-2)])
	lpsi <- as.double(param[(npar-1):npar])
	omega<-as.double(param[npar+1])
  	theta <- work <- as.double(rep(0,n))
  	g.beta <- d.beta <- d.beta1<-d.beta2 <- der <- as.double(rep(0,npar-2))
  	g.lpsi1<-g.lpsi2 <- as.double(0)
  	d.beta <- matrix(as.double(0),3,npar-2)
	gvar<-as.double(0)
 	x <- matrix(as.double(X),nrow=n, ncol=npar-2)
  	der <- as.double(rep(0,npar-2))
  	db <- matrix(as.double(0),3,npar-2)
  	db1<- matrix(as.double(0),4,npar-2)
  	db2<- matrix(as.double(0),5,npar-2)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

  	result <- .Fortran("gint",g.beta,g.lpsi1,g.lpsi2,gvar,
                  x,theta,work,y,lpsi,beta,bt,d.beta,d.beta1,d.beta2,
                  der,db,db1,db2,omega, npar,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")

	return(c(result[[1]],result[[2]],result[[3]],result[[4]]))}

	loglik2 <- function(param, X, y,integrate)
	{
	npar <- as.integer(length(param)-1)
	beta <- as.double(param[1:(npar-2)])
	bt <- as.double(param[1:(npar-2)])
	lpsi <- as.double(param[(npar-1):npar])
	omega<-as.double(param[npar+1])
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
	x <-matrix(as.double(X),nrow=n,ncol=npar-2)
	theta <- work <-as.double(rep(0,n))
	logL <- as.double(0)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

 		results <- .Fortran("integ",logL,bt,beta,lpsi,omega,npar,
		x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")


	return(results[[1]])}
	
	nparam <- as.integer(length(parameters)-1)
	omega1<-parameters[nparam+1]
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	dbeta<-as.double(rep(0,nparam-2))
	dlog.psi1<-0
	dlog.psi2<-0
	dvar<-0
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik2(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
	grad<-gradient2(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)

	for (j in 1:(nparam-2))
	{
	dbeta[j]<-dbeta[j]+counts[i]*(grad[j]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
	}
	dlog.psi1<-dlog.psi1+counts[i]*(grad[nparam-1]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
	dlog.psi2<-dlog.psi2+counts[i]*(grad[nparam]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
 	#using the chain rule
	dvar<-dvar+counts[i]*(grad[nparam+1]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))*exp(omega1)
	k1<-k2+1
	}
gr<-c(dbeta,dlog.psi1,dlog.psi2,dvar)
return(-gr)}

# ******************* MAIN PROGRAM *******************************
#
	call <- match.call()
#	vect.time <- F
	if(missing(data) || !is.data.frame(data))
		stop("a data.frame must be supplied")
	if(is.null(names(data)))
		stop("objects in data.frame must have a name")
	expr1 <- terms(formula, data=data)
	expr <- attr(expr1, "variables")
	var.names <- all.vars(expr)
	response <- all.vars(expr)[1]

	if(any(is.na(match(var.names, names(data)))))
		stop("Variables in formula not contained in the data.frame")

if(!missing(time)) {Time<-as.vector(data[[time]])}
if (missing(time)) { if (all (is.na(match(names(data), "time")))) stop ("time must be defined")
			else  Time<-as.vector(data$time)}
				
if(!missing(id)) {id<-as.vector(data[[id]])}
if (missing(id)){ if (all(is.na(match(names(data), "id")))) stop ("id must be defined")
		else id<-as.vector(data$id)}

# select subset if necessary
	if(!missing(subSET)) {id1 <- eval(substitute(subSET), data)
			      data<-subset(data, id1)}

	if(!missing(aggregate)) {f.name <- deparse(substitute(aggregate))
  			f.value <- as.factor(data[[f.name]])}
	

	if(missing(aggregate)) {f.value<-as.factor(0)}

#returns data of a subset
subset.data<-data

	ti.repl<-as.vector(0)
	i1<-1
	i2<-1
	for (i in 1:(length(data[[response]])-1))
	{ 
		if (id[i]==id[i+1])
		{ 	i2<-i2+1
			ti.repl[i1]<-i2}
		else {	ti.repl[i1]<-i2
			i1<-i1+1
			i2<-1}
	}

	n.cases <- length(ti.repl)
	n.tot<-cumsum(ti.repl)[n.cases]
	n.time<-length(unique(Time))
	ni.cases <- length(ti.repl)
	pos.ind<-cumsum(ti.repl)

	counts<-as.vector(0)

	if(is.na(match("counts", names(data))))
			counts <- data$counts <- rep(1, n.cases)
		else   {for (i in 1:n.cases)
			{counts[i]<-data$counts[pos.ind[i]]}
		}

	final.data <- data
	var.names <- c(var.names, "counts")
	data <- data[var.names]
	n.var <- length(data)
	Y.resp <- as.vector(data[[response]])
	if((all(Y.resp[!is.na(Y.resp)] == 1 | Y.resp[!is.na(Y.resp)] == 0)) == FALSE)
		stop("Unfeasible values of response variable: must be 0,1,NA")

# ********** creation of individual profile according to NA patterns *******************

	data2<-data
		
	if (dependence=="MC2"|| dependence=="MC2R")

	{final.data <- na.discrete.replace2(frame=data,  n.times=n.time, ti.repl=ti.repl)}

	else if (dependence=="ind"|| dependence=="indR"||dependence=="MC1"|| dependence=="MC1R")

	{final.data <- na.discrete.replace1(frame=data,  n.times=n.time, ti.repl=ti.repl)}


if (length(final.data$bad.ind)>=1)
	{data<-final.data$data
	data<-data[-final.data$badlines,]
	data2<-data2[-final.data$badlines,]
	counts<-counts[-final.data$bad.ind]
	ti.repl<-ti.repl[-final.data$bad.ind]
	n.cases <- length(ti.repl)
	n.tot<-cumsum(ti.repl)[n.cases]}
else {data<-final.data$data}


# ********** design matrices creation *******************
	# define a plausible starting point for the optimizer if not given
		data1 <- na.omit(data2)
		data1.resp <- data1[, response]
		data1.resp <- log((data1.resp + 0.5)/(1.5 - data1.resp))
		data1[, c(response)] <- data1.resp

		if (dependence=="MC1")	 init<-0
		else  if (dependence=="MC1R")  init<-c(0,0)
		else  if (dependence=="MC2")    init<-c(0,0)
		else  if (dependence=="MC2R")   init<-c(0,0,0)
		else  if (dependence=="indR")  init<-0

		if(is.null(start) && dependence!="ind")
			start <- c(lm(formula, data1, weights = counts)$coefficients, init)
		else if(!is.null(start) && dependence!="ind")
			start <- c(lm(formula, data1, weights = counts)$coefficients, start)
		else if (dependence=="ind") start <- c(lm(formula, data1, weights = counts)$coefficients)

if (any(is.na(start))) stop("starting values produced by lm contains NA")

	id.not.na<-rep(TRUE,n.tot)
	X <- model.matrix(expr1, data, contrasts)
	names.output <- dimnames(X)[[2]]
	sum.ti <- sum(ti.repl)
	data <- list(ti.repl, data[[response]], counts)
	data2<-list(ti.repl, data2[[response]], counts)
	p <- dim(X)[2] + 1
	pr<-prob<-F.aux<-as.double(rep(0,length(data[[2]])))

	if (dependence=="ind")
	{	if(trace)	cat("\t log.likelihood\n")
	temp <-optim(par= start, fn = logL.bin0,  gr= gradlogL.bin0, method=method, 
	data = data, X = X, trace=trace,control=control)}

	else  if (dependence=="indR")
	{	if(trace)	cat("\n omega\t log.likelihood\n")
	temp <- optim(par = start, fn =logL.bin0Ifm,gr = gradlogL.bin0Ifm,  method=method, 
	data = data, X = X, integrate=integrate, trace=trace,control=control)}

	else if (dependence=="MC1")
	{	if(trace)	cat("\n log.psi1\t log.likelihood\n")
	temp <-optim(par= start, fn = logL.bin1mf ,  gr= gradlogL.bin1mf, method=method, 
	data = data, X = X, trace=trace,control=control)}
	else  if (dependence=="MC1R")
	{	if(trace)	cat("\n log.psi1\t omega\t log.likelihood\n")
	temp <- optim(par = start, fn =logL.bin1Ifm,gr = gradlogL.binI1fm,  method=method, 
	data = data, X = X, integrate=integrate, trace=trace,control=control)}
	else  if (dependence=="MC2")
	{	if(trace)	cat("\n log.psi1\t log.psi2\t log.likelihood\n")
	temp <- optim(par= start, fn =logL.bin2mf, gr =gradlogL.bin2mf, method=method,
	 data = data, X = X, trace = trace,control=control)}
	else  if (dependence=="MC2R")
	{	if(trace)	cat("\n log.psi1\t log.psi2\t omega\t log.likelihood\n")
	temp <- optim(par = start, fn= logL.bin2Ifm, gr = gradlogL.binI2fm,method=method, 
	data = data, X = X,integrate=integrate, trace = trace,control=control)}

	coefficients <- temp$par
	log.lik <-  - temp$value
	if (trace) 
	cat("Convergence reached. Computing the information matrix now\n")

	 if (dependence=="ind")
	Info <- num.info(coefficients, "gradlogL.bin0", X, data)
	else  if (dependence=="indR")
	Info <- num.infoI(coefficients, "gradlogL.bin0Ifm", X, data, integrate=integrate)
	else if (dependence=="MC1")
	Info <- num.info(coefficients, "gradlogL.bin1mf", X, data)
	else  if (dependence=="MC1R")
	Info <- num.infoI(coefficients, "gradlogL.binI1fm", X, data, integrate=integrate)
	else  if (dependence=="MC2")
	Info <- num.info(coefficients, "gradlogL.bin2mf", X, data)
	else  if (dependence=="MC2R")
	Info <- num.infoI(coefficients, "gradlogL.binI2fm", X, data, integrate=integrate)

	se <- matrix(sqrt(diag(solve(Info))), ncol = 1)
	coefficients <- matrix(coefficients, ncol = 1)
	if (dependence=="ind")
	dimnames(coefficients) <- dimnames(se) <- list(names.output, " ")
	else  if (dependence=="indR")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "omega"), " ")
	else if (dependence=="MC1")
	dimnames(coefficients) <- dimnames(se) <- list(c(names.output, "log.psi1"), " ")
	else  if (dependence=="MC1R")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "log.psi1","omega"), " ")
	else if (dependence=="MC2")
	dimnames(coefficients) <- dimnames(se) <- list(c(names.output, "log.psi1","log.psi2"), " ")
	else if (dependence=="MC2R")
	dimnames(coefficients) <- dimnames(se) <- list(c(names.output, "log.psi1","log.psi2","omega" ), " ")

	covariance <- solve(Info)
	cr<- diag(1/sqrt(diag(covariance)))
	correlation <- cr %*% covariance %*% cr

	if (dependence=="ind")
	{dimnames(covariance) <- list(names.output, names.output)
	dimnames(correlation) <- list(names.output, names.output)}

	else  if (dependence=="indR")
	{dimnames(covariance) <- list(c(names.output, "omega"), c(names.output, "omega"))
	dimnames(correlation) <- list(c(names.output, "omega"), c(names.output, "omega"))}

	else if (dependence=="MC1")
	{dimnames(covariance) <- list(c(names.output, "log.psi1"), c(names.output, "log.psi1"))
	dimnames(correlation) <- list(c(names.output, "log.psi1"), c(names.output, "log.psi1"))}
	else  if (dependence=="MC1R")
	{dimnames(covariance) <- list(c(names.output, "log.psi1","omega"), c(names.output, 	"log.psi1","omega"))
	dimnames(correlation) <- list(c(names.output, "log.psi1","omega"), c(names.output, 	"log.psi1","omega"))}
	else if (dependence=="MC2")
	{dimnames(covariance) <- list(c(names.output, "log.psi1","log.psi2"), c(names.output,"log.psi1","log.psi2"))
	dimnames(correlation) <- list(c(names.output, "log.psi1","log.psi2"), c(names.output,"log.psi1","log.psi2"))}
	else if (dependence=="MC2R")
	{dimnames(covariance) <- list(c(names.output, "log.psi1","log.psi2","omega"), c(names.output, 		"log.psi1","log.psi2","omega"))
	dimnames(correlation) <- list(c(names.output, "log.psi1","log.psi2","omega"), c(names.output, 		"log.psi1","log.psi2","omega"))}

#### To compute estimated transition probabilities
 	if (dependence=="ind")
	{pr <- logL.bin0.aux (parameters=coefficients, X=X, data=data, trace=trace)
	prob<-pr$pij}
	else  if (dependence=="indR")
	{pr <- logL.bin0I.aux (parameters=coefficients, X=X, data=data2, trace=trace)
	prob<-pr$pij}
	else if (dependence=="MC1")
	{pr <- logL.bin1mf.aux (parameters=coefficients, X=X, data=data, trace=trace)
	prob<-pr$pij}
	else  if (dependence=="MC1R")
	{pr <- logL.bin1.aux (parameters=coefficients, X=X, data=data2, trace=trace)
	prob<-pr$pij}
	else  if (dependence=="MC2")
	{pr <- logL.bin2mf.aux(parameters=coefficients, X=X, data=data, trace=trace)
	prob<-pr$pij}
	else  if (dependence=="MC2R")
	{pr <- logL.bin2.aux(parameters=coefficients, X=X, data=data2, trace=trace)
	prob<-pr$pij}

#### To compute fitted values 
	Fitted <- rep(NA, n.tot)
 	if (dependence=="ind"|dependence=="MC1"|dependence=="MC2")
	{Fitted[id.not.na] <- X %*% coefficients[1:(p - 1)]}
	else  if (dependence=="indR")
	{F.aux<- logL.bin0I.aux (parameters=coefficients, X=X, data=data2, trace=trace)
	Fitted<-F.aux$fit}
	else  if (dependence=="MC1R")
	{F.aux<- logL.bin1.aux (parameters=coefficients, X=X, data=data2, trace=trace)
	Fitted<-F.aux$fit}
	else  if (dependence=="MC2R")
	{F.aux <- logL.bin2.aux(parameters=coefficients, X=X, data=data2, trace=trace)
	Fitted<-F.aux$fit}

	ncoef<-length(coefficients)
	aic<-(2*temp$value+2*ncoef)
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data2[[2]]
	counts<-data[[3]]
	residuals<-wts<-freq<-as.double(rep(0,length(y)))
	ncounts<-sum(counts[1:n.cases])
	residuals<-(y-prob)/sqrt(prob*(1-prob))

	Fitted <- plogis(Fitted)
	Fitted[is.na(y)] <- NA
	res<-as.double(rep(0,n.time))
	
	for (j in 1: n.time)
	{ k2<-0
	soma.n<-soma.d<-0
     	for ( i in 1:n.cases)
	{ 
	k3<-k2+j
 	if(!is.na(y[k3]))
     		{
      		soma.n<-soma.n+(y[k3]-prob[k3])
  		soma.d<-soma.d+(prob[k3]*(1-prob[k3]))
       		}
	k2<-cumti.repl[i]
  	}
	res[j]<-soma.n/sqrt(soma.d)
	}

	y.matrix<-matrix(y,ncol=n.time,byrow=TRUE)
	y.av<-apply(y.matrix,2,mean,na.rm=TRUE)
	Fitted.matrix<-matrix(Fitted,ncol=n.time,byrow=TRUE)
	Fitted.av<-apply(Fitted.matrix,2,mean,na.rm=TRUE)
	prob.matrix<-matrix(prob,ncol=n.time,byrow=TRUE)


bl<- new("bild", coefficients = coefficients, se = se, covariance =covariance, correlation=correlation, 
	log.likelihood=- temp$value, message = temp$convergence, n.cases=n.cases, ni.cases=ni.cases, aic=aic,    
      	residuals=residuals, s.residuals=res, ind.probability=prob, prob.matrix=prob.matrix, Fitted=Fitted, 
	Fitted.av=Fitted.av, Time=Time, model.matrix=X, y.matrix=y.matrix, subset.data=subset.data,
	y.av=y.av, f.value=f.value, call=call)

}



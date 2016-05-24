
setClass("summary.cold", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric",call="language"))

setClass("cold", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric",  
		 Fitted="numeric", Fitted.av="numeric", Time="numeric", model.matrix= "matrix", y.matrix="matrix",
		subset.data="data.frame",final.data="data.frame", y.av="numeric", f.value="factor", data.id="numeric",call="language"))


setGeneric("getAIC",def=function(object) standardGeneric("getAIC"))


setGeneric("getLogLik",def=function(object) standardGeneric("getLogLik"))


cold<-function(formula = formula(data), data, time,id, subSET, aggregate=FALSE, start = NULL, trace = FALSE,
dependence="ind", method="BFGS", control=coldControl(), 
integrate=coldIntegrate())
{
# *****************DEFINITION OF INTERNAL FUNCTIONS ******************

na.discrete.replace <- function(frame,  n.times, ti.repl)
	{
	vars <- names(frame)
	names(vars) <- vars
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	for(j in 1:length(vars)) 
	{k1<-1
	for (i in 1:n.cases)
	{k2<-cumti.repl[i]
	x <- frame[[j]][k1:k2]
	pos <- is.na(x)
	if(any(pos))
			if(j == 1) x[pos] <- -1
			else x[pos] <- x[1]
			
	frame[[j]][k1:k2]<-x
	k1<-k2+1
	}
	}
		
#print(c("na.discrete.replace", frame))
	return(data=frame)
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

#################################
#					  #
# função log_likelihood		  #
#					  #
#################################

logL.pss0 <- function(parameters, X, data, trace)
{

loglik1<- function(param, X, y, trace)
{#calculate logLi(beta/bi) for each individual

npar <- as.integer(length(param)+1)
beta <- as.double(param[1:(npar-1)])

y[is.na(y)]<-(-1)
y <- as.integer(y)
n <- as.integer(length(y)) 
x <-matrix(as.double(X),nrow=n,ncol=npar-1)
theta <- work <- as.double(rep(0,n))
logL <- L<-as.double(0)


		link <- as.integer(1)
		m <- max(y)
		fact <- rep(1, m + 1)
		if(m > 0)
			{for(i in 2:(m + 1))
			fact[i] <- fact[i - 1] * (i - 1)}
			fact <- as.double(fact)

		results <- .Fortran("psslik0",logL,beta,
			npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
				

	return(results[[1]])
}


ti.repl<-data[[1]]
cumti.repl<-cumsum(ti.repl)
n.cases<- length(ti.repl)
y<-data[[2]]
counts<-data[[3]]
logL<-0
k1<-1

	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)

	logL<-logL+counts[i]*z
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")

return(-logL)}


# ------------------------------------------------------------------- 

LogL.pss1 <- function(parameters, X, data, trace)
{

loglik1<- function(param, X, y, trace)
{#calculate logLi(beta/bi) for each individual

npar <- as.integer(length(param))
beta <- as.double(param[1:(npar-1)])
rho<- as.double(param[npar])
y[is.na(y)]<-(-1)
y <- as.integer(y)
n <- as.integer(length(y)) 
x <-matrix(as.double(X),nrow=n,ncol=npar-1)
theta <- work <-as.double(rep(0,n))
logL <- L<-as.double(0)



		link <- as.integer(1)
		m <- max(y)
		fact <- rep(1, m + 1)
		if(m > 0)
			{for(i in 2:(m + 1))
			fact[i] <- fact[i - 1] * (i - 1)}
			
		fact <- as.double(fact)

		results <- .Fortran("psslik",logL,beta,rho,
			npar,x,y,theta,work,n,fact,link,PACKAGE="cold")

	return(results[[1]])
	}
		if(trace)	cat(paste(format(parameters[length(parameters)], digit=4), collapse=" "), "\t")

ti.repl<-data[[1]]
cumti.repl<-cumsum(ti.repl)
n.cases<- length(ti.repl)
y<-data[[2]]
counts<-data[[3]]
logL<-0
k1<-1

npar <- as.integer(length(parameters))
rho<- as.double(parameters[npar])

if (rho < 0 |  rho >1 ) 
{ logL<-NaN
		if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
return(NaN)}

else if (rho > 0 &  rho < 1 ) 
{
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)

	logL<-logL+counts[i]*z
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
return(-logL)}
}

######################################################################## começa novo
# ------------------------------------------------------------------- 
	
LogL.pss0I<- function(parameters, X, data, integrate, trace)
{

loglik1 <- function(param, X, y,integrate, trace)
{
npar <-as.integer(length(param))
beta<- as.double(param[1:(npar-1)])
bt<- as.double(param[1:(npar-1)])
omega<-as.double(param[npar])
y[is.na(y)]<-(-1)
y<- as.integer(y)
n <- as.integer(length(y)) 
x<-matrix(as.double(X),nrow=n,ncol=npar-1)
theta<- work<- as.double(rep(0,n))
logL <- as.double(0)

li<-as.double(integrate$li)
ls<-as.double(integrate$ls)
epsabs<-as.double(integrate$epsabs)
epsrel<-as.double(integrate$epsrel)
limit<-as.integer(integrate$limit)
key<-as.integer(integrate$key)

m <- max(y)
	


	if(m >20)
	{	li<-as.double(-0.01)
		ls<-as.double(0.01)
		key<-as.integer(4)
	}

	link <- as.integer(1)

	 results <- .Fortran("intp0",logL,bt,beta,omega,
	npar,link,m,x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="cold")

		logL<-results[[1]]

return(results[[1]])}
	
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

nparam<-length(parameters)
omega1<-as.double(parameters[nparam])
ti.repl<-data[[1]]
cumti.repl<-cumsum(ti.repl)
n.cases<- as.integer(length(ti.repl))
y<-data[[2]]
counts<-data[[3]]
logL1<-0
k1<-1

	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate,trace=trace)
	#logL1 gives the log-likelihood
	logL1<-logL1+counts[i]*log(z*(1/(sqrt(2*pi)*exp(omega1/2))))
	k1<-k2+1
	}
		if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
return(-logL1)}


############### to compute fitted

LogL.pss0I.aux<- function(parameters, X, data, trace)
{
	loglik1 <- function(param, X, y,trace)
	{
	
	npar <-as.integer(length(param))
	beta<- as.double(param[1:(npar-1)])
	bt<- as.double(param[1:(npar-1)])
	omega<-as.double(param[npar])
	y<- as.integer(y)
	n <- as.integer(length(y)) 
	x<-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta<- work<-as.double(rep(0,n))
	logL <- as.double(0)

	eta<-i.fit<- fit<-as.vector(npar-1)
	eta<-x%*%beta
	m<-glm(as.numeric(y)~offset(eta), family=poisson)
	bi<-coef(m)
	beta[1]<-beta[1]+bi
	y[is.na(y)]<-(-1)
	y<- as.integer(y)


		link <- as.integer(1)
		m <- max(y)
		fact <- rep(1, m + 1)
		if(m > 0) 
			{for(i in 2:(m + 1))
			fact[i] <- fact[i - 1] * (i - 1)}
			
		fact <- as.double(fact)

	results <- .Fortran("psslik0",logL,beta,
			npar,x,y,theta,work,n,fact,link,PACKAGE="cold")


	i.fit<-x%*%beta

	return(list(loglik=results[[1]], fit=i.fit))}
	
		if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=4), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

	nparam<-length(parameters)
	omega1<-as.double(parameters[nparam])
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	fitted<-as.double(rep(0,length(y)))
	logL1<-as.double(0)
	k1<-1

	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],trace=trace)
	fitted[k1:k2]<-z$fit
	k1<-k2+1
	}

return(fit=fitted)}


# ------------------------------------------------------------------- acaba novo


LogL.pss1I<- function(parameters, X, data, integrate, trace)
{

loglik1 <- function(param, X, y,integrate, trace)
{
npar <-as.integer(length(param)-1)
beta<- as.double(param[1:(npar-1)])
bt<- as.double(param[1:(npar-1)])
rho<-as.double(param[npar])
omega<-as.double(param[npar+1])
y[is.na(y)]<-(-1)
y<- as.integer(y)
n <- as.integer(length(y)) 
x<-matrix(as.double(X),nrow=n,ncol=npar-1)
theta<- work<- as.double(rep(0,n))
logL <- as.double(0)

li<-as.double(integrate$li)
ls<-as.double(integrate$ls)
epsabs<-as.double(integrate$epsabs)
epsrel<-as.double(integrate$epsrel)
limit<-as.integer(integrate$limit)
key<-as.integer(integrate$key)

m <- max(y)

	link <- as.integer(1)

	 results <- .Fortran("intp",logL,bt,beta,rho,omega,
	npar,link,m,x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="cold")

		logL<-results[[1]]

return(results[[1]])}
	
		if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=4), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

nparam<-length(parameters)
omega1<-as.double(parameters[nparam])
ti.repl<-data[[1]]
cumti.repl<-cumsum(ti.repl)
n.cases<- as.integer(length(ti.repl))
y<-data[[2]]
counts<-data[[3]]
logL1<-0
k1<-1
npar <- as.integer(length(parameters))
rho<- as.double(parameters[npar-1])

if (rho < 0 |  rho >1 ) 
{ logL<-NaN
		if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
return(NaN)}

else if (rho > 0 &  rho < 1 ) 
{
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate,trace=trace)
	#logL1 gives the log-likelihood
	logL1<-logL1+counts[i]*log(z*(1/(sqrt(2*pi)*exp(omega1/2))))
	k1<-k2+1
	}

		if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
return(-logL1)}
}

############### to compute fitted

LogL.pss1I.aux<- function(parameters, X, data, trace)
{
	loglik1 <- function(param, X, y,trace)
	{
	
	npar <-as.integer(length(param)-1)
	beta<- as.double(param[1:(npar-1)])
	bt<- as.double(param[1:(npar-1)])
	rho<-as.double(param[npar])
	omega<-as.double(param[npar+1])
	y<- as.integer(y)
	n <- as.integer(length(y)) 
	x<-matrix(as.double(X),nrow=n,ncol=npar-1)
	theta<- work<-as.double(rep(0,n))
	logL <- as.double(0)

	eta<-i.fit<- fit<-as.vector(npar-1)
	eta<-x%*%beta
	m<-glm(as.numeric(y)~offset(eta), family=poisson)
	bi<-coef(m)
	beta[1]<-beta[1]+bi
	y[is.na(y)]<-(-1)
	y<- as.integer(y)

		link <- as.integer(1)
		m <- max(y)
		fact <- rep(1, m + 1)
		if(m > 0)
			{for(i in 2:(m + 1))
			fact[i] <- fact[i - 1] * (i - 1)}
			fact <- as.double(fact)
	
	results <- .Fortran("psslik",logL,beta,rho,
			npar,x,y,theta,work,n,fact,link,PACKAGE="cold")


	i.fit<-x%*%beta

	return(list(loglik=results[[1]], fit=i.fit))}
	
		if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=4), collapse=" "), "\t")
		if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")

	nparam<-length(parameters)
	omega1<-as.double(parameters[nparam])
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data[[2]]
	counts<-data[[3]]
	fitted<-as.double(rep(0,length(y)))
	logL1<-as.double(0)
	k1<-1
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],trace=trace)
	fitted[k1:k2]<-z$fit
	k1<-k2+1
	}

return(fit=fitted)}

# ------------------------------------------------------------------- 


#################################
#				#
# função gradientlog_likelihood #
#				#
#################################


gradlogL.pss0<- function(parameters, X,data, trace)
{

gradient <- function(param, X, y)
		{
		npar <- as.integer(length(param)+1)
		beta <- as.double(param[1:(npar-1)])
		rho <- as.double(0)
		y[is.na(y)]<-(-1)
		y <- as.integer(y)
		n <- as.integer(length(y)) 
		theta <- work <- as.double(rep(0, n))
		grad <- as.double(rep(0, npar))
		x <-matrix(as.double(X),nrow=n,ncol=npar-1)
		m <- max(y)
		fact <- rep(1, m + 1)
		if(m > 0)
			{for(i in 2:m + 1)
			fact[i] <- fact[i - 1] * (i - 1)}
			fact <- as.double(fact)
			link <- as.integer(1)
		
		result <- .Fortran("pssgrd0",grad,beta,rho,
			npar,x,y,theta,work,n,fact,link,PACKAGE="cold")

	for (i in 1:length(grad))
	{if  (result[[1]][i]=="NaN" ) result[[1]][i]<-0}

		return(result[[1]])}

	ti.repl <- data[[1]]
	cumti.repl <- cumsum(ti.repl)
	n.cases <- length(ti.repl)
	y <- data[[2]]
	counts <- data[[3]]
	gr<-0
	k1 <- 1
		
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	u<-gradient (param=parameters,X=X[k1:k2,], y=y[k1:k2])

	gr<-gr+counts[i] * u
	k1<-k2+1
	}

	grad.f<-gr[1:length(parameters)]


 return(-grad.f)}

####---------------------------------------------------------

gradLogL.pss1<- function(parameters, X,data, trace)
{

gradient <- function(param,  X, y)
		{
		npar <- as.integer(length(param))
		beta <- as.double(param[1:(npar-1)])
		rho <- as.double(param[npar])
		y[is.na(y)] <- (-1)
		y <- as.integer(y)
		n <- as.integer(length(y)) 
		theta <- work <- as.double(rep(0, n))
		grad <- as.double(rep(0, npar))
		x <-matrix(as.double(X),nrow=n,ncol=npar-1)
		m <- max(y)
		fact <- rep(1, m + 1)
		if(m > 0)
			{for(i in 2:m + 1)
			fact[i] <- fact[i - 1] * (i - 1)}
			fact <- as.double(fact)
			link <- as.integer(1)
		
		result <- .Fortran("pssgrd",grad,beta,rho,
			npar,x,y,theta,work,n,fact,link,PACKAGE="cold")

	for (i in 1:length(grad))
	{if  (result[[1]][i]=="NaN" ) result[[1]][i]<-0}

		return(result[[1]])}

	ti.repl <- data[[1]]
	cumti.repl <- cumsum(ti.repl)
	n.cases <- length(ti.repl)
	y <- data[[2]]
	counts <- data[[3]]
	gr<-0
	k1 <- 1
		
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	u<-gradient (param=parameters,X=X[k1:k2,], y=y[k1:k2])
	gr<-gr+counts[i] * u
	k1<-k2+1
	}

	 return(-gr)}

####-------------------------------------------

gradLogL.pss0I <- function(parameters, X,data,integrate,trace)
{

gradient1 <-  function(param,X,y,integrate)
	{
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
  	npar <- as.integer(length(param))
	beta <- as.double(param[1:(npar-1)])
	bt <- as.double(param[1:(npar-1)])
	omega<-as.double(param[npar])
  	theta <- work <- as.double(rep(0,n))
  	grad<- as.double(rep(0,npar-1))
	gvar<-as.double(0)
 	x <- matrix(as.double(X),nrow=n, ncol=npar-1)
	m <- max(y)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

	if(m >20)
	{	li<-as.double(-0.000001)
		ls<-as.double(0.000001)
		key<-as.integer(3)
	}

		link <- as.integer(1)
		
	result <- .Fortran("gintp0",grad,gvar,bt,beta,omega,npar,link,
     		m,x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="cold")
	

	return(c(result[[1]],result[[2]]))	}


loglik1 <- function(param, X, y,integrate)
{
npar <-as.integer(length(param))
beta<- as.double(param[1:(npar-1)])
bt<- as.double(param[1:(npar-1)])
omega<-as.double(param[npar])
y[is.na(y)]<-(-1)
y<- as.integer(y)
n <- as.integer(length(y)) 
x<-matrix(as.double(X),nrow=n,ncol=npar-1)
theta<- work<- as.double(rep(0,n))
logL <- as.double(0)
m <- max(y)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

	if(m >20)
	{	li<-as.double(-0.0001)
		ls<-as.double(0.0001)
		key<-as.integer(3)
	}


link <- as.integer(1)

 	results <- .Fortran("intp0",logL,bt,beta,omega,npar,link,m,
		x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="cold")

return(results[[1]])}

	
nparam <- as.integer(length(parameters))
omega1<-parameters[nparam]
ti.repl<-data[[1]]
cumti.repl<-cumsum(ti.repl)
n.cases<- length(ti.repl)
y<-data[[2]]
counts<-data[[3]]
dgr<-as.double(rep(0,nparam-1))
dvar<-0
k1<-1
	

	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
	
	num<-gradient1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)


	for (j in 1:(nparam-1))
	{
	if (is.na(num[j]) | is.na(num[j]-Inf)) num[j]<-0
	dgr[j]<-dgr[j]+counts[i]*(num[j]/z)	
	}

 #using the chain rule
	if (is.na(num[nparam]) | is.na(num[nparam]-Inf)) num[nparam]<-0
	dvar<-dvar+counts[i]*(num[nparam]/z)*exp(omega1)
	k1<-k2+1
	}

gr<-c(dgr,dvar)


return(-gr)}

####-------------------------------------------

gradLogL.pss1I <- function(parameters, X,data,integrate,trace)
{

gradient1 <-  function(param,X,y,integrate)
	{
	y[is.na(y)]<-(-1)
	y <- as.integer(y)
	n <- as.integer(length(y)) 
  	npar <- as.integer(length(param)-1)
	beta <- as.double(param[1:(npar-1)])
	bt <- as.double(param[1:(npar-1)])
	rho <- as.double(param[npar])
	omega<-as.double(param[npar+1])
  	theta <- work <- as.double(rep(0,n))
  	grad<- as.double(rep(0,npar))
	gvar<-as.double(0)
 	x <- matrix(as.double(X),nrow=n, ncol=npar-1)
	m <- max(y)


	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

	if(m >20)
	{	li<-as.double(-0.000001)
		ls<-as.double(0.000001)
		key<-as.integer(3)
	}

		link <- as.integer(1)
		
	result <- .Fortran("gintp",grad,gvar,bt,beta,rho,omega,npar,link,
     		m,x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="cold")
	


	return(c(result[[1]],result[[2]]))
	}

loglik1 <- function(param, X, y,integrate)
{
npar <-as.integer(length(param)-1)
beta<- as.double(param[1:(npar-1)])
bt<- as.double(param[1:(npar-1)])
rho<-as.double(param[npar])
omega<-as.double(param[npar+1])
y[is.na(y)]<-(-1)
y<- as.integer(y)
n <- as.integer(length(y)) 
x<-matrix(as.double(X),nrow=n,ncol=npar-1)
theta<- work<- as.double(rep(0,n))
logL <- as.double(0)
m <- max(y)

	li<-as.double(integrate$lig)
	ls<-as.double(integrate$lsg)
	epsabs<-as.double(integrate$epsabs)
	epsrel<-as.double(integrate$epsrel)
	limit<-as.integer(integrate$limit)
	key<-as.integer(integrate$key)

	if(m >20)
	{	li<-as.double(-0.0001)
		ls<-as.double(0.0001)
		key<-as.integer(3)
	}

	link <- as.integer(1)

 	results <- .Fortran("intp",logL,bt,beta,rho,omega,npar,link,m,
		x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="cold")

return(results[[1]])}

	
nparam <- as.integer(length(parameters)-1)
omega1<-parameters[nparam+1]
ti.repl<-data[[1]]
cumti.repl<-cumsum(ti.repl)
n.cases<- length(ti.repl)
y<-data[[2]]
counts<-data[[3]]
dgr<-as.double(rep(0,nparam))
dvar<-0
k1<-1
	
	for (i in 1:n.cases)
	{
	k2<-cumti.repl[i]
	
	z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
	
	num<-gradient1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)

	for (j in 1:(nparam))
	{
	if (is.na(num[j]) | is.na(num[j]-Inf)) num[j]<-0
	dgr[j]<-dgr[j]+counts[i]*(num[j]/z)	
	}

 #using the chain rule
	if (is.na(num[nparam+1]) | is.na(num[nparam+1]-Inf)) num[nparam+1]<-0
	dvar<-dvar+counts[i]*(num[nparam+1]/z)*exp(omega1)

	k1<-k2+1
	}
gr<-c(dgr,dvar)
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
	y1 <- Y.resp[!is.na(Y.resp)]
	if((all(y1 >= 0) | all(y1[y1 != 0] == as.integer(y1[y1 != 0]))) == FALSE)
	stop("Unfeasible values of response variable: must be non negative integers or NA"	)


# ********** creation of individual profile according to NA patterns *******************

#### tenho dúvidas pq foi adaptado do bild (para caso de missing data) e n tenho certeza se está correto!

	data2<-data
	final.data <- na.discrete.replace(frame=data,  n.times=n.time, ti.repl=ti.repl)
	
	data<-final.data

# ********** design matrices creation *******************
	# define a plausible starting point for the optimizer if not given
		data1 <- na.omit(data2)
		data1.resp <- data1[, response]
#		data1.resp <- log(data1.resp + 1)
		data1[, c(response)] <- data1.resp

		if (dependence=="AR1")	 init<-0.5
		else  if (dependence=="AR1R")  init<-c(0.5,0)
		else  if (dependence=="indR")  init<-0

		if(is.null(start) && dependence!="ind")
			start <- c(glm(formula, data1,family=poisson, weights=counts)$coefficients, init)
		else if(!is.null(start) && dependence!="ind")
			start <- c(glm(formula, data1,family=poisson, weights=counts)$coefficients, start)
		else if (dependence=="ind") start <- c(glm(formula, data1, family=poisson, weights=counts)$coefficients)

if (any(is.na(start))) stop("starting values produced by glm contains NA")

	id.not.na<-rep(TRUE,n.tot)
	X <- model.matrix(expr1, data, contrasts)
	names.output <- dimnames(X)[[2]]
	sum.ti <- sum(ti.repl)
	data <- list(ti.repl, data[[response]], counts)
	data2<-list(ti.repl, data2[[response]], counts)
	p <- dim(X)[2] + 1
	F.aux<-as.double(rep(0,length(data[[2]])))


	if (dependence=="ind")
	{	if(trace)	cat("\t log.likelihood\n")
	temp<-optim(par= start, fn=logL.pss0, gr=gradlogL.pss0, method=method, 
	data = data, X = X, trace=trace,control=control)}
	else  if (dependence=="indR")
	{	if(trace)	cat("\n omega\t log.likelihood\n")
	temp <- optim(par = start, fn =LogL.pss0I,gr = gradLogL.pss0I,  method=method, 
	data = data, X = X, integrate=integrate, trace=trace,control=control)}
	else if (dependence=="AR1")
	{	if(trace)	cat("\n rho\t log.likelihood\n")
	temp <-optim(par= start, fn = LogL.pss1,  gr=gradLogL.pss1, method=method, 
	data = data, X = X, trace=trace,control=control)}
	else  if (dependence=="AR1R")
	{	if(trace)	cat("\n rho\t omega\t log.likelihood\n")
	temp <- optim(par = start, fn =LogL.pss1I,gr = gradLogL.pss1I,  method=method, 
	data = data, X = X, integrate=integrate, trace=trace,control=control)}

	coefficients <- temp$par
	log.lik <-  - temp$value
	if (trace) 
	cat("Convergence reached. Computing the information matrix now\n")

	if (dependence=="ind")
	Info <- num.info(coefficients, "gradlogL.pss0", X, data)
	else  if (dependence=="indR")
	Info <- num.infoI(coefficients, "gradLogL.pss0I", X, data, integrate=integrate) 
	else if (dependence=="AR1")
	Info <- num.info(coefficients, "gradLogL.pss1", X, data)
	else  if (dependence=="AR1R")
	Info <- num.infoI(coefficients, "gradLogL.pss1I", X, data, integrate=integrate)

	se <- matrix(sqrt(diag(solve(Info))), ncol = 1)
	coefficients <- matrix(coefficients, ncol = 1)
	if (dependence=="ind")
	dimnames(coefficients) <- dimnames(se) <- list(names.output, " ")
	else  if (dependence=="indR")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "omega"), " ")
	else if (dependence=="AR1")
	dimnames(coefficients) <- dimnames(se) <- list(c(names.output, "rho"), " ")
	else  if (dependence=="AR1R")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "rho","omega"), " ")
	
	covariance <- solve(Info)
	cr<- diag(1/sqrt(diag(covariance)))
	correlation <- cr %*% covariance %*% cr


	if (dependence=="ind")
	{dimnames(covariance) <- list(names.output, names.output)
	dimnames(correlation) <- list(names.output, names.output)}
	else  if (dependence=="indR")
	{dimnames(covariance) <- list(c(names.output, "omega"), c(names.output, "omega"))
	dimnames(correlation) <- list(c(names.output, "omega"), c(names.output, "omega"))}
	else if (dependence=="AR1")
	{dimnames(covariance) <- list(c(names.output, "rho"), c(names.output, "rho"))
	dimnames(correlation) <- list(c(names.output, "rho"), c(names.output, "rho"))}
	else  if (dependence=="AR1R")
	{dimnames(covariance) <- list(c(names.output, "rho","omega"), c(names.output, 	"rho","omega"))
	dimnames(correlation) <- list(c(names.output, "rho","omega"), c(names.output, 	"rho","omega"))}


#### To compute fitted values 
	Fitted <- rep(NA, n.tot)
 	if (dependence=="ind"|dependence=="AR1")
	{Fitted[id.not.na] <- X %*% coefficients[1:(p - 1)]}
	else  if (dependence=="indR")
	{Fitted<- LogL.pss0I.aux (parameters=coefficients, X=X, data=data2, trace=trace)}
	else  if (dependence=="AR1R")
	{Fitted<- LogL.pss1I.aux (parameters=coefficients, X=X, data=data2, trace=trace)}

	ncoef<-length(coefficients)
	aic<-(2*temp$value+2*ncoef)
	y<-data2[[2]]

	Fitted <- exp(Fitted)
	Fitted[is.na(y)] <- NA

	y.matrix<-matrix(y,ncol=n.time,byrow=TRUE)
	y.av<-apply(y.matrix,2,mean,na.rm=TRUE)
	Fitted.matrix<-matrix(Fitted,ncol=n.time,byrow=TRUE)
	Fitted.av<-apply(Fitted.matrix,2,mean,na.rm=TRUE)


cl<- new("cold", coefficients = coefficients, se = se, covariance =covariance, correlation=correlation, 
	log.likelihood=- temp$value, message = temp$convergence, n.cases=n.cases, ni.cases=ni.cases, aic=aic,    
      	Fitted=Fitted, Fitted.av=Fitted.av, Time=Time, model.matrix=X, 
	y.matrix=y.matrix, subset.data=subset.data, final.data=final.data,y.av=y.av, f.value=f.value, 
	data.id=id,call=call)

}



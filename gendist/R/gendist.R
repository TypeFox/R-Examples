####################################################################
# Mixture Model
####################################################################
dmixt	<- function(x, phi, spec1, arg1, spec2, arg2, log=FALSE){
	if(phi<0) {stop("phi must be of non-negative values")}
	f1	<- function(d,z) {do.call(paste(d, spec1, sep=""), c(list(z), arg1) )}
	f2	<- function(d,z) {do.call(paste(d, spec2, sep=""), c(list(z), arg2) )}
	pdf	<- (1/(1+phi))*f1("d", x) + (1-(1/(1+phi)))*f2("d",x)
	tt 	<- if(log) {log(pdf)} else {pdf}
	return(tt)
}

pmixt	<- function(q, phi, spec1, arg1, spec2, arg2, lower.tail=TRUE, log.p=FALSE){
	if(phi<0) {stop("phi must be of non-negative values")}
	f1	<- function(d,z) {do.call(paste(d, spec1, sep=""), c(list(z), arg1) )}
	f2	<- function(d,z) {do.call(paste(d, spec2, sep=""), c(list(z), arg2) )}
	cdf	<- (1/(1+phi))*f1("p", q) + (1-(1/(1+phi)))*f2("p",q)
	tt	<- if(lower.tail) {cdf} else {1-cdf}
	tt 	<- if(log.p) {log(tt)} else {tt}
	return(tt)
}

qmixt	<- function(p, phi, spec1, arg1, spec2, arg2, interval=c(0,100), lower.tail=TRUE, log.p=FALSE){
	if(phi<0) {stop("phi must be of non-negative values")}
	p	<- if(lower.tail) {p} else {1-p}
	qf	<- function(p) {uniroot(function(q) {ifelse(p <=0, 0,
		   ifelse(log.p,
		   pmixt(q, phi, spec1, arg1, spec2, arg2, lower.tail, log.p)-log(p),
		   pmixt(q, phi, spec1, arg1, spec2, arg2, lower.tail, log.p)-p))}, interval,tol=1e-10)$root}
	qf	<- Vectorize(qf)
	tt	<- qf(p)
	return(tt)
}

rmixt	<- function(n, phi, spec1, arg1, spec2, arg2, interval=c(0,100)){
	if(phi<0) {stop("phi must be of non-negative values")}
	rf	<- qmixt(runif(n), phi, spec1, arg1, spec2, arg2, interval)
	return(rf)
}


####################################################################
# Composite Model
####################################################################
dcomposite<-function(x, spec1, arg1, spec2, arg2, initial=1, log=FALSE)
{	
	fun1		<- function(d,z) {do.call(paste(d, spec1, sep=""), c(list(z), arg1) )}
	fun2		<- function(d,z) {do.call(paste(d, spec2, sep=""), c(list(z), arg2) )}
	initial	<- initial

	ratio		<- function(p){-log(fun1("d",p)/fun2("d",p))}
	theta		<- nlm(ratio,initial)$estimate
	phi		<- ( fun1("d",theta)*(1-fun2("p",theta) ) / ( fun2("d",theta)*fun1("p",theta) ) )

	pdf		<- ifelse( x<=theta,
					(1/(1+phi)) / fun1("p",theta)*fun1("d",x),
					(phi/(1+phi)) / (1-fun2("p",theta))*fun2("d",x)
			   )
	tt 		<- if(log) {log(pdf)} else {pdf}
	return(tt)
}

################################################
pcomposite<-function(q, spec1, arg1, spec2, arg2, initial=1, lower.tail=TRUE, log.p=FALSE)
{	
	fun1		<- function(d,z) {do.call(paste(d, spec1, sep=""), c(list(z), arg1) )}
	fun2		<- function(d,z) {do.call(paste(d, spec2, sep=""), c(list(z), arg2) )}

	ratio		<- function(p){-log(fun1("d",p)/fun2("d",p))}
	theta		<- nlm(ratio,initial)$estimate
	phi			<- ( fun1("d",theta)*(1-fun2("p",theta) ) / ( fun2("d",theta)*fun1("p",theta) ) )

	cdf			<- ifelse( q<=theta,
					(1/(1+phi)) / fun1("p",theta)*fun1("p",q),
					(1/(1+phi)) + (phi/(1+phi)) / (1-fun2("p",theta))*(fun2("p",q)-fun2("p",theta))
				)
	tt	<- if(lower.tail) {cdf} else {1-cdf}
	tt 	<- if(log.p) {log(tt)} else {tt}
	return(tt)
}

################################################
qcomposite<-function(p, spec1, arg1, spec2, arg2, initial=1, lower.tail=TRUE, log.p=FALSE)
{	
	p		<- if(lower.tail) {p} else {1-p}
	fun1		<- function(d,z) {do.call(paste(d, spec1, sep=""), c(list(z), arg1) )}
	fun2		<- function(d,z) {do.call(paste(d, spec2, sep=""), c(list(z), arg2) )}

	ratio		<- function(p){-log(fun1("d",p)/fun2("d",p))}
	theta		<- nlm(ratio,initial)$estimate
	phi		<- ( fun1("d",theta)*(1-fun2("p",theta) ) / ( fun2("d",theta)*fun1("p",theta) ) )

	
	qf		<- function(p) {ifelse( p<=0, 0,
				ifelse(log.p,
					ifelse( p<=(1/(1+phi)),
						fun1("q", log(p)*(1+phi)*fun1("p",theta)),
						fun2("q", fun2("p",theta) + (log(p)*(1+phi)-1)*(1-fun2("p",theta))/phi)
					),
					ifelse( p<=(1/(1+phi)),
						fun1("q", p*(1+phi)*fun1("p",theta)),
						fun2("q", fun2("p",theta) + (p*(1+phi)-1)*(1-fun2("p",theta))/phi)
					)
				)
			   )}
	qf		<- Vectorize(qf)
	tt		<- qf(p)
	return(tt)
}

################################################
rcomposite		<- function(n, spec1, arg1, spec2, arg2, initial=1){
	rf		<- qcomposite(runif(n), spec1, arg1, spec2, arg2, initial)
	return(rf)
}


####################################################################
# Folded Model
####################################################################
####################################################
dfolded <- function(x, spec, arg, log=FALSE)
{	
	fun	<- function(d,z) {do.call(paste(d, spec, sep=""), c(list(z), arg) )}	
	pdf	<- ifelse(x<0,0,fun("d",x)+fun("d",-x))
	tt	<- if(log) {log(pdf)} else {pdf}
	return(tt)
}

####################################################
pfolded <- function(q, spec, arg, lower.tail=TRUE, log.p=FALSE)
{	
	fun	<- function(d,z) {do.call(paste(d, spec, sep=""), c(list(z), arg) )}	
	cdf	<- ifelse(q<=0,0,fun("p", q)-fun("p", -q))
	
	tt	<- if(lower.tail) {cdf} else {1-cdf}
	tt 	<- if(log.p) {log(tt)} else {tt}
	return(tt)
}

####################################################
qfolded <- function(p, spec, arg, interval=c(0,100), lower.tail=TRUE, log.p=FALSE)
{	
	p	<- if(lower.tail) {p} else {1-p}
	qf	<- function(p) {uniroot(function(q) {ifelse(p <=0, 0,
		   	ifelse(log.p,
		  		pfolded(q, spec, arg, lower.tail, log.p)-log(p),
		   		pfolded(q, spec, arg, lower.tail, log.p)-p))}, interval,tol=1e-10)$root}
	qf	<- Vectorize(qf)
	tt	<- qf(p)
	return(tt)
}

####################################################
rfolded <- function(n, spec, arg, interval=c(0,100))
{
	rf	<- qfolded(runif(n),spec, arg, interval)
	return(rf)
}


#############################
#	Skewed Models
#############################
dskew <- function(x, spec1, arg1, spec2, arg2, log=FALSE)
{	
	fun1	<- function(d,z) {do.call(paste(d, spec1 , sep=""), c(list(z), arg1) )}
	fun2	<- function(d,z) {do.call(paste(d, spec2 , sep=""), c(list(z), arg2) )}

	pdf	<- 2*fun1("d",x)*fun2("p",x)
	tt	<- if(log) {log(pdf)} else {pdf}
	return(tt)
}

####################################################
pskew <- function(q, spec1, arg1, spec2, arg2, lower.tail=TRUE, log.p=FALSE)
{	
	fun1	<- function(d,z) {do.call(paste(d, spec1 , sep=""), c(list(z), arg1) )}
	fun2	<- function(d,z) {do.call(paste(d, spec2 , sep=""), c(list(z), arg2) )}

	cdf	<- function(q) integrate(function(q) {2*fun1("d", q)*fun2("p", q)}, lower=-Inf, upper=q)$value
	cdf	<- Vectorize(cdf)
	cdf	<- cdf(q)

	tt	<- if(lower.tail) {cdf} else {1-cdf}
	tt 	<- if(log.p) {log(tt)} else {tt}
	return(tt)
}

####################################################
qskew <- function(p, spec1, arg1, spec2, arg2, interval=c(1,10), lower.tail=TRUE, log.p=FALSE)
{	
	p	<- if(lower.tail) {p} else {1-p}
	fun1	<- function(d,z) {do.call(paste(d, spec1 , sep=""), c(list(z), arg1) )}
	fun2	<- function(d,z) {do.call(paste(d, spec2 , sep=""), c(list(z), arg2) )}

	qf	<- function(p) {uniroot(function(q) {ifelse(p <=0, 0,
		   ifelse(log.p,
		   pskew(q, spec1, arg1, spec2, arg2, lower.tail, log.p)-log(p),
		   pskew(q, spec1, arg1, spec2, arg2, lower.tail, log.p)-p))}, interval,tol=1e-10)$root}
	qf	<- Vectorize(qf)
	tt	<- qf(p)
	return(tt)
}

####################################################
rskew <- function(n, spec1, arg1, spec2, arg2, interval=c(1, 10))
{
	rf	<- qskew(runif(n), spec1, arg1, spec2, arg2, interval)
	return(rf)
}


#############################
#	Arc Tan Models
#############################
darctan <- function(x, alpha, spec, arg, log=FALSE)
{	
	if(alpha<=0) {stop("alpha must be of positive values")}
	fun	<- function(d,z) {do.call(paste(d, spec, sep=""), c(list(z), arg) )}	
	pdf	<- (1/atan(alpha)) * (alpha*fun("d",x))/(1+ (alpha*(1-fun("p",x)))^2)
	tt	<- if(log) {log(pdf)} else {pdf}
	return(tt)
}

####################################################
parctan <- function(q, alpha, spec, arg, lower.tail=TRUE, log.p=FALSE)
{	
	if(alpha<=0) {stop("phi must be of non-negative values")}
	fun	<- function(d,z) {do.call(paste(d, spec, sep=""), c(list(z), arg) )}	
	cdf	<- 1- atan(alpha*(1-fun("p", q))) / atan(alpha)
	
	tt	<- if(lower.tail) {cdf} else {1-cdf}
	tt 	<- if(log.p) {log(tt)} else {tt}
	return(tt)
}

####################################################
qarctan <- function(p, alpha, spec, arg, lower.tail=TRUE, log.p=FALSE)
{
	if(alpha<=0) {stop("alpha must be of positive values")}
	p	<- if(lower.tail) {p} else {1-p}
	fun	<- function(d,z) {do.call(paste(d, spec, sep=""), c(list(z), arg) )}		

	qf	<- function(p) {ifelse(p <=0, 0,
		   	ifelse(log.p,
		  		fun("q", 1 - (tan(atan(alpha)*(1-log(p))))/alpha ),
		   		fun("q", 1 - (tan(atan(alpha)*(1-p)))/alpha )
				) )}
	qf	<- Vectorize(qf)
	tt	<- qf(p)
	return(tt)
}

####################################################
rarctan <- function(n, alpha, spec, arg)
{
	if(alpha<=0) {stop("alpha must be of positive values")}
	rf	<- qarctan(runif(n), alpha, spec, arg)
	return(rf)
}



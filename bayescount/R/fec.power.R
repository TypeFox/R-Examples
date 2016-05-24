fec.power.limits <- function(meanepg=200, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, coeffvarrep=0.4, coeffvarind=0.3, coeffvargroup=0.7, true.sample=FALSE, lower.limit=NA, upper.limit=NA, iterations=100000, power = 0.95, confidence = 0.99, feedback=FALSE, forcesim=FALSE){
	
# Power is the required power, confidence is the confidence in this power when using simulations
if(power >= 1) stop("Required power must be < 1")

if(confidence >= 1) stop("Confidence must be < 1")
conf <- (1-confidence)/2
lci <- 0+conf
uci <- 1-conf

if(replicates < 1 | animals < 1) stop("Specified values for animals and replicates must be greater than 0")
if(!is.na(lower.limit) & !is.na(upper.limit)) stop("One or both of lower.limit or upper.limit must be non-fixed (NA)")

fix.lower <- !is.na(lower.limit)
fix.upper <- !is.na(upper.limit)

fixed.lower <- lower.limit
fixed.upper <- upper.limit

target <- power

if(feedback){
if(any(.Platform$GUI == c("AQUA", "Rgui"))){
#feedback <- FALSE
warning("Printing the progress of the function using the GUI versions of R may massively increase the time taken, I suggest setting feedback=FALSE or using a command line version of R instead")
}
##### CHECK ON WINDOWS / UNIX
}

# If we have 1 animal and true.sample=T then we can forget about cv.group and use nbpower
if(animals==1 & true.sample==TRUE) coeffvargroup <- 10^-10
if(coeffvargroup > (coeffvarrep+coeffvarind)/10) approximate <- FALSE else approximate <- TRUE
if(forcesim) approximate <- FALSE

if(animals==1 & true.sample==FALSE) warning("NOTE:  The power calculated is for the population from which the animal is derived, to calculate the power for the true mean of this individual use true.sample=TRUE")

lowerl <- lower.limit
upperl <- upper.limit

start <- Sys.time()

if(!approximate){
	if(true.sample){
		
		out <- .C("poweranalysissamplefixed", as.numeric(meanepg), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(coeffvarrep), as.numeric(coeffvarind), as.numeric(coeffvargroup), as.integer(iterations), as.integer(feedback), numeric(iterations), PACKAGE="bayescount")

		lo <- length(out)
		meancounts <- out[[lo]]
		
	}else{
		out <- .C("poweranalysispopulationfixed", as.numeric(meanepg), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(coeffvarrep), as.numeric(coeffvarind), as.numeric(coeffvargroup), as.integer(iterations), as.integer(feedback), numeric(iterations), PACKAGE="bayescount")

		lo <- length(out)
		meancounts <- out[[lo]]
		
	}
	
	mcs <- meancounts
	
	f <- function(limit){
		prob <- limit
		#limits <- HPDinterval(mcmc(mcs), prob=prob)[1,]
		limits <- quantile(mcs, probs=c(0+((1-prob)/2), 1-((1-prob)/2)))
		nin <- sum(mcs <= limits[2] & mcs >= limits[1])
		nout <- length(mcs)-nin
		med <- qbeta(0.5, nin+1, nout+1)
		return(med)
	}

	if(fix.lower){
	f <- function(limit){
		#limits <- HPDinterval(mcmc(mcs), prob=prob)[1,]
		nin <- sum(mcs <= limit & mcs >= fixed.lower)
		nout <- length(mcs)-nin
		med <- qbeta(0.5, nin+1, nout+1)
		return(med)
	}
	}
	if(fix.upper){
	f <- function(limit){
		#limits <- HPDinterval(mcmc(mcs), prob=prob)[1,]
		nin <- sum(mcs <= fixed.upper & mcs >= -limit)
		nout <- length(mcs)-nin
		med <- qbeta(0.5, nin+1, nout+1)
		return(med)
	}
	}

	
	limits <- c(0,1)
	if(fix.upper) limits <- c(-Inf,0)
	if(fix.lower) limits <- c(0,Inf)
	bsres <- binary.search(f, target, limits)
	if(bsres$status!="OK"){
		if(bsres$status!="Absolute value not possible but this is closest") stop(bsres$status)
	}

	if(fix.upper) limits <- c(-bsres$value, fixed.upper)
	if(fix.lower) limits <- c(fixed.lower, bsres$value)
	if(!fix.upper & !fix.lower) limits <- quantile(mcs, prob=c(0+((1-bsres$value)/2), 1-((1-bsres$value)/2)))

	nin <- sum(mcs <= limits[2] & mcs >= limits[1])
	nout <- length(mcs)-nin
	
	power <- qbeta(c(lci,0.5,uci), nin+1, nout+1)
	
	names(limits) <- c("lower.limit", "upper.limit")
	names(power) <- c(paste("lower.", confidence*100, "%.estimate", sep=""), "median.estimate", paste("upper.", confidence*100, "%.estimate", sep=""))
	
	return(list(limits=limits, power=power))
}else{
	#Embed within wrapper function (that calculates cvw from rep, sample, g.faeces etc and works out if approximation to multiple animals can be used):

	if(true.sample==FALSE | animals > 1){
		# If we're using this with true.sample=FALSE then cvgroup must be so small compared to cvind and cvrep that we can just pretend that we're drawing all replicates from the same animal (with a slightly larger cv and animalmean=populationmean):
		coeffvarind <- sqrt(coeffvarind^2 + coeffvargroup^2 + coeffvarind^2*coeffvargroup^2)
		replicates <- replicates*animals
	}
	
	coeff.var <- sqrt(coeffvarind^2 + (coeffvarrep^2)/g.faeces + (coeffvarind^2 * (coeffvarrep^2)/g.faeces))
	eggs.counted <- replicates*meanepg*sensitivity
	shape <- replicates / coeff.var^2
	
	f <- function(limit){
		accuracy <- limit
		
		upper.tol <- replicates*sensitivity*(meanepg+meanepg*accuracy)
		lower.tol <- replicates*sensitivity*(meanepg-meanepg*accuracy)
		
		return(pnbinom(upper.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) - pnbinom(lower.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) + suppressWarnings(dnbinom(lower.tol, size=shape, mu=eggs.counted)))
		
	}	
	if(fix.lower){
	f <- function(limit){
		upper.tol <- replicates*sensitivity*(limit)
		lower.tol <- replicates*sensitivity*(fixed.lower)
		
		return(pnbinom(upper.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) - pnbinom(lower.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) + suppressWarnings(dnbinom(lower.tol, size=shape, mu=eggs.counted)))
	}
	}
	if(fix.upper){
	f <- function(limit){
		upper.tol <- replicates*sensitivity*(fixed.upper)
		lower.tol <- replicates*sensitivity*(-limit)
		
		return(pnbinom(upper.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) - pnbinom(lower.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) + suppressWarnings(dnbinom(lower.tol, size=shape, mu=eggs.counted)))
	}
	}

	
	limits <- c(0,Inf)
	if(fix.upper) limits <- c(-Inf,0)
	if(fix.lower) limits <- c(0,Inf)
	bsres <- binary.search(f, target, limits)
	if(bsres$status!="OK"){
		if(bsres$status!="Absolute value not possible but this is closest") stop(bsres$status)
	}

	if(fix.upper) limits <- c(-bsres$value, fixed.upper)
	if(fix.lower) limits <- c(fixed.lower, bsres$value)
	if(!fix.upper & !fix.lower) limits <- c((meanepg-meanepg*bsres$value), (meanepg+meanepg*bsres$value))
	
	power <- replicate(3, bsres$objective)
	
	names(limits) <- c("lower.limit", "upper.limit")
	names(power) <- c("absolute.value", "absolute.value", "absolute.value")
	
	return(list(limits=limits, power=power))
	
}


}



fec.power <- function(meanepg=200, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, coeffvarrep=0.4, coeffvarind=0.3, coeffvargroup=0.7, true.sample=FALSE, accuracy=0.1, lower.limit=meanepg*(1-accuracy), upper.limit=meanepg*(1+accuracy), maxiterations=1000000, precision=2, confidence = 0.99, feedback=FALSE, forcesim=FALSE){

if(confidence >= 1) stop("Confidence must be < 1")
conf <- (1-confidence)/2
lci <- 0+conf
uci <- 1-conf

if(feedback){
if(any(.Platform$GUI == c("AQUA", "Rgui"))){
#feedback <- FALSE
warning("Printing the progress of the function using the GUI versions of R may massively increase the time taken, I suggest setting feedback=FALSE or using a command line version of R instead")
}
##### CHECK ON WINDOWS / UNIX
}

if(replicates < 1 | animals < 1) stop("Specified values for animals and replicates must be greater than 0")

# If we have 1 animal and true.sample=T then we can forget about cv.group and use nbpower
if(animals==1 & true.sample==TRUE) coeffvargroup <- 10^-10
if(coeffvargroup > (coeffvarrep+coeffvarind)/10) approximate <- FALSE else approximate <- TRUE
if(forcesim) approximate <- FALSE

if(animals==1 & true.sample==FALSE) warning("NOTE:  The power calculated is for the population from which the animal is derived, to calculate the power for the true mean of this individual use true.sample=TRUE")

lowerl <- lower.limit
upperl <- upper.limit

start <- Sys.time()

if(approximate){
	
		if(true.sample==FALSE){
			# If we're using this with true.sample=FALSE then cvgroup must be so small compared to cvind and cvrep that we can just pretend that we're drawing all replicates from the same animal (with a slightly larger cv and animalmean=populationmean):
			coeffvarind <- sqrt(coeffvarind^2 + coeffvargroup^2 + coeffvarind^2*coeffvargroup^2)
			replicates <- replicates*animals
		}
		
		# Produces the same results as poweranalysispopulation (only for animal=1 and/or cv.group<<cv.other)
		
		coeff.var <- sqrt(coeffvarind^2 + (coeffvarrep^2)/g.faeces + (coeffvarind^2 * (coeffvarrep^2)/g.faeces))

		eggs.counted <- replicates*meanepg*sensitivity

		upper.tol <- replicates*sensitivity*(upperl)
		lower.tol <- replicates*sensitivity*(lowerl)

		shape <- replicates / coeff.var^2
		
		nbpower <- pnbinom(upper.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) - pnbinom(lower.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) + suppressWarnings(dnbinom(lower.tol, size=shape, mu=eggs.counted))
		
		time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
		
		output = list(roundedci=replicate(3, round(nbpower, precision)), ci=replicate(3, nbpower))
	
}else{
	if(true.sample){
		if(is.na(precision)){
			out <- .C("poweranalysissamplefixed", as.numeric(meanepg), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(coeffvarrep), as.numeric(coeffvarind), as.numeric(coeffvargroup), as.integer(maxiterations), as.integer(feedback), numeric(maxiterations), PACKAGE="bayescount")

			lo <- length(out)
			meancounts <- out[[lo]]
			nin <- sum(meancounts >= lowerl & meancounts <= upperl)
			nout <- length(meancounts)-nin
			time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
			
			output=list(roundedci=round(qbeta(c(lci, 0.5, uci), nin+1, nout+1), precision), ci=qbeta(c(lci, 0.5, uci), nin+1, nout+1), within=nin, without=nout, total=nin+nout)
			
		}else{
			out <- .C("poweranalysissample", as.numeric(meanepg), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(coeffvarrep), as.numeric(coeffvarind), as.numeric(coeffvargroup), as.numeric(lowerl), as.numeric(upperl), as.integer(maxiterations), as.integer(precision), as.numeric(lci), as.numeric(uci), as.integer(feedback), as.integer(0), as.integer(0), PACKAGE="bayescount")

			lo <- length(out)
			time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
			
			output=list(roundedci = round(qbeta(c(lci, 0.5, uci), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), precision), ci=qbeta(c(lci, 0.5, uci), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), within=out[[lo-1]], without=out[[lo]]-out[[lo-1]], total=out[[lo]])
		}
	
	}else{
		
		# coeffvarrep is the cv between subsamples taken from the same replicate and is effected by gfaeces (calculation done in C)
		# coeffvarind is cv between replicates and is effected by replicates.  This actually means it is EITHER the combined cv between samples and days within an animal OR just the cv between samples if taken from a single faecal pile (in which case the cv between days is actually part of the cv between animals and will consequently make things a bit wrong due to the true sample effect).  Doesn't reflect subsamples (what I usually call replicates) as this is exactly the same as changing the sensitivity and g.faeces.
		# coeffvargroup is cv between animals and is always calculated explicitly since there is no shortcut.  If we have only 1 animal or coeffvargroup=0, then powernb works perfectly.

		if(is.na(precision)){
			out <- .C("poweranalysispopulationfixed", as.numeric(meanepg), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(coeffvarrep), as.numeric(coeffvarind), as.numeric(coeffvargroup), as.integer(maxiterations), as.integer(feedback), numeric(maxiterations), PACKAGE="bayescount")
			
			lo <- length(out)
			meancounts <- out[[lo]]
			nin <- sum(meancounts >= lowerl & meancounts <= upperl)
			nout <- length(meancounts)-nin
			time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
			
			output=list(roundedci=round(qbeta(c(lci, 0.5, uci), nin+1, nout+1), precision), ci=qbeta(c(lci, 0.5, uci), nin+1, nout+1), within=nin, without=nout, total=nin+nout)
			
			
		}else{
			out <- .C("poweranalysispopulation", as.numeric(meanepg), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(coeffvarrep), as.numeric(coeffvarind), as.numeric(coeffvargroup), as.numeric(lowerl), as.numeric(upperl), as.integer(maxiterations), as.integer(precision), as.numeric(lci), as.numeric(uci), as.integer(feedback), as.integer(0), as.integer(0), PACKAGE="bayescount")

			lo <- length(out)
			time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
			
			output=list(roundedci = round(qbeta(c(lci, 0.5, uci), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), precision), ci=qbeta(c(lci, 0.5, uci), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), within=out[[lo-1]], without=out[[lo]]-out[[lo-1]], total=out[[lo]])
		}
	}
}

if(approximate){
	names(output$roundedci) <- c("absolute.value", "absolute.value", "absolute.value")
	names(output$ci) <- c("absolute.value", "absolute.value", "absolute.value")
}else{
	names(output$roundedci) <- c(paste("lower.", confidence*100, "%", sep=""), "median", paste("upper.", confidence*100, "%", sep=""))
	names(output$ci) <- c(paste("lower.", confidence*100, "%", sep=""), "median", paste("upper.", confidence*100, "%", sep=""))
}

return(output)
}

FEC.power.limits <- fec.power.limits
FEC.precision <- fec.power.limits
fec.precision <- fec.power.limits
count.precision <- fec.precision
FEC.power <- fec.power
count.power <- fec.power

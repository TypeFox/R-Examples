fecrt.power.limits <- function(meanepg=200, reduction = 95, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, pre.coeffvarrep=0.4, pre.coeffvarind=0.3, pre.coeffvargroup=0.7, post.coeffvarrep=0.4, post.coeffvarind=0.3, post.coeffvargroup=0.7, true.sample=FALSE, lower.limit=NA, upper.limit=NA, iterations=100000, power = 0.95, confidence = 0.99, feedback=FALSE){
	
# Power is the required power, confidence is the confidence in this power when using simulations
if(power >= 1) stop("Required power must be < 1")

if(confidence >= 1) stop("Confidence must be < 1")
conf <- (1-confidence)/2
lci <- 0+conf
uci <- 1-conf

if(replicates < 1 | animals < 1) stop("Specified values for animals and replicates must be greater than 0")

if(!is.na(lower.limit) & !is.na(upper.limit)) stop("One or both of lower.limit or upper.limit must be non-fixed (NA)")

fix.upper <- !is.na(lower.limit)
fix.lower <- !is.na(upper.limit)

fixed.lower = lowerl <- (100-upper.limit)/100
fixed.upper = upperl <- (100-lower.limit)/100
delta <- (100-reduction)/100

tlowerl <- lowerl
tupperl <- upperl
tdelta <- delta
if(is.na(tlowerl)) tlowerl <- 0.5
if(is.na(tupperl)) tupperl <- 0.5
if(is.na(tdelta)) tdelta <- 0.5

if(tlowerl < 0 | tupperl < 0 | tdelta < 0 | tlowerl > 1 | tupperl > 1 | tdelta > 1) stop("reduction, lower.limit and upper.limit (if supplied) must all be between 0 and 100(%)")

target <- power

if(feedback){
if(any(.Platform$GUI == c("AQUA", "Rgui"))){
#feedback <- FALSE
warning("Printing the progress of the function using the GUI versions of R may massively increase the time taken, I suggest setting feedback=FALSE or using a command line version of R instead")
}
##### CHECK ON WINDOWS / UNIX
}


if(animals==1 & true.sample==FALSE) warning("NOTE:  The power calculated is for the population from which the animal is derived, to calculate the power for the true mean of this individual use true.sample=TRUE")

lowerl <- lower.limit
upperl <- upper.limit

delta <- (100-reduction)/100

if(delta < 0 | delta > 1) stop("Supplied value for 'reduction' must be between 0 and 100(%)")

start <- Sys.time()

if(true.sample){
	
	out <- .C("fecrtpowersamplefixed", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.integer(iterations), as.integer(feedback), numeric(iterations), PACKAGE="bayescount")

	lo <- length(out)
	meanreds <- out[[lo]]
	
}else{
	out <- .C("fecrtpowerpopulationfixed", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.integer(iterations), as.integer(feedback), numeric(iterations), PACKAGE="bayescount")

	lo <- length(out)
	meanreds <- out[[lo]]
	
}

mcs <- meanreds

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
if(fix.upper) limits <- c(-1,0)
if(fix.lower) limits <- c(0,1)

bsres <- binary.search(f, target, limits)
if(bsres$status!="OK"){
	if(bsres$status!="Absolute value not possible but this is closest"){
		limits <- c(0,1)
		if(fix.upper) limits[2] <- fixed.upper
		if(fix.lower) limits[1] <- fixed.lower
		nin <- sum(mcs <= limits[1] & mcs >= limits[2])
		nout <- length(mcs)-nin
		med <- round(qbeta(0.5, nin+1, nout+1), 2)
		stop(paste("It is not possible to achieve the desired power with any tolerance; maximum power =", med))
	}
}

if(fix.upper) limits <- c(-bsres$value, fixed.upper)
if(fix.lower) limits <- c(fixed.lower, bsres$value)
if(!fix.upper & !fix.lower) limits <- quantile(mcs, prob=c(0+((1-bsres$value)/2), 1-((1-bsres$value)/2)))

nin <- sum(mcs <= limits[2] & mcs >= limits[1])
nout <- length(mcs)-nin

#limits <- paste(prettyround((100*(1-limits))[2:1], 2), "%", sep="")
limits <- (100*(1-limits))[2:1]

power <- qbeta(c(lci,0.5,uci), nin+1, nout+1)

names(limits) <- c("lower.limit", "upper.limit")
names(power) <- c(paste("lower.", confidence*100, "%.estimate", sep=""), "median.estimate", paste("upper.", confidence*100, "%.estimate", sep=""))

return(list(limits=limits, power=power))

}




fecrt.power <- function(meanepg=200, reduction = 80, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, pre.coeffvarrep=0.4, pre.coeffvarind=0.3, pre.coeffvargroup=0.7, post.coeffvarrep=0.4, post.coeffvarind=0.3, post.coeffvargroup=0.7, true.sample=FALSE, lower.limit = 0, upper.limit=90, maxiterations=1000000, precision=2, confidence = 0.99, feedback=FALSE){

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

if(animals==1 & true.sample==FALSE) warning("NOTE:  The power calculated is for the population from which the animal is derived, to calculate the power for the true mean of this individual use true.sample=TRUE")

lowerl <- (100-upper.limit)/100
upperl <- (100-lower.limit)/100
delta <- (100-reduction)/100

if(lowerl < 0 | upperl < 0 | delta < 0 | lowerl > 1 | upperl > 1 | delta > 1) stop("reduction, lower.limit and upper.limit must all be between 0 and 100(%)")

start <- Sys.time()
if(true.sample){
	if(is.na(precision)){
		out <- .C("fecrtpowersamplefixed", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.integer(maxiterations), as.integer(feedback), numeric(maxiterations), PACKAGE="bayescount")

		lo <- length(out)
		meancounts <- out[[lo]]
				
		nin <- sum(meancounts >= lowerl & meancounts <= upperl)
		nout <- length(meancounts)-nin
		time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
		
		output=list(roundedci=round(qbeta(c(lci, 0.5, uci), nin+1, nout+1), precision), ci=qbeta(c(lci, 0.5, uci), nin+1, nout+1), within=nin, without=nout, total=nin+nout)
		
	}else{
		out <- .C("fecrtpowersample", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.numeric(lowerl), as.numeric(upperl), as.integer(maxiterations), as.integer(precision), as.numeric(lci), as.numeric(uci), as.integer(feedback), as.integer(0), as.integer(0), PACKAGE="bayescount")

		lo <- length(out)
		time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
		
		output=list(roundedci = round(qbeta(c(lci, 0.5, uci), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), precision), ci=qbeta(c(lci, 0.5, uci), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), within=out[[lo-1]], without=out[[lo]]-out[[lo-1]], total=out[[lo]])
	}

}else{
	
	# coeffvarrep is the cv between subsamples taken from the same replicate and is effected by gfaeces (calculation done in C)
	# coeffvarind is cv between replicates and is effected by replicates.  This actually means it is EITHER the combined cv between samples and days within an animal OR just the cv between samples if taken from a single faecal pile (in which case the cv between days is actually part of the cv between animals and will consequently make things a bit wrong due to the true sample effect).  Doesn't reflect subsamples (what I usually call replicates) as this is exactly the same as changing the sensitivity and g.faeces.
	# coeffvargroup is cv between animals and is always calculated explicitly since there is no shortcut.  If we have only 1 animal or coeffvargroup=0, then powernb works perfectly.

	if(is.na(precision)){
		out <- .C("fecrtpowerpopulationfixed", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.integer(maxiterations), as.integer(feedback), numeric(maxiterations), PACKAGE="bayescount")
		
		lo <- length(out)
		meancounts <- out[[lo]]
				
		nin <- sum(meancounts >= lowerl & meancounts <= upperl)
		nout <- length(meancounts)-nin
		time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
		
		output=list(roundedci=round(qbeta(c(lci, 0.5, uci), nin+1, nout+1), precision), ci=qbeta(c(lci, 0.5, uci), nin+1, nout+1), within=nin, without=nout, total=nin+nout)
		
		
	}else{
		out <- .C("fecrtpowerpopulation", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.numeric(lowerl), as.numeric(upperl), as.integer(maxiterations), as.integer(precision), as.numeric(lci), as.numeric(uci), as.integer(feedback), as.integer(0), as.integer(0), PACKAGE="bayescount")

		lo <- length(out)
		time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
		
		output=list(roundedci = round(qbeta(c(lci, 0.5, uci), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), precision), ci=qbeta(c(lci, 0.5, uci), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), within=out[[lo-1]], without=out[[lo]]-out[[lo-1]], total=out[[lo]])
	}
}

names(output$roundedci) <- c(paste("lower.", confidence*100, "%", sep=""), "median", paste("upper.", confidence*100, "%", sep=""))
names(output$ci) <- c(paste("lower.", confidence*100, "%", sep=""), "median", paste("upper.", confidence*100, "%", sep=""))

return(output)
}

FECRT.power.limits <- fecrt.power.limits
FECRT.power <- fecrt.power
fecrt.precision <- fecrt.power.limits
FECRT.precision <- fecrt.power.limits

bayescount.single <- function(data = stop("Data must be specified"), model="ZILP", alt.prior = FALSE, adjust.zi.mean = FALSE, raw.output = FALSE, likelihood=FALSE, ...){


passthrough <- list(...)
if(is.null(passthrough$max.time)) passthrough$max.time <- "1hr"
if(is.null(passthrough$interactive)) passthrough$interactive <- FALSE
if(is.null(passthrough$plots)) passthrough$plots <- FALSE

arguments <- formals(runjags::autorun.jags)
#newargs <- formals(runjags::add.summary)
#arguments <- c(arguments, newargs[!names(newargs) %in% names(arguments)])
arguments <- arguments[! names(arguments) %in% c('...', 'runjags.object')]

for(i in 1:length(passthrough)){
	arguments[[names(passthrough)[i]]] <- passthrough[[i]]
}

jags <- eval(arguments$jags)
silent.jags <- eval(arguments$silent.jags)

test <- testjags(jags, silent=TRUE)
if(test[[2]][1]==FALSE){
	cat("Unable to call JAGS using '", jags, "' - try specifying the path to the JAGS binary as the jags argument\n", sep="")
	stop("Unable to call JAGS")
}


if(alt.prior!=FALSE & (any(model==c("IP", "SP", "ZISP")))){
	alt.prior <- FALSE
	#  Warning printed in run.model
}

if(class(raw.output)=="list"){
	record.name <- raw.output$name
	record.setname <- raw.output$setname
	record.chains <- TRUE
	raw.output <- FALSE
}else{
	record.chains <- FALSE
}


adjust.mean <- adjust.zi.mean

model <- toupper(model)

model <- switch(model, P="SP", ZIP="ZISP", model)

true.model <- model
if(any(model==paste("G", 1:10, sep=""))){
	true.model <- model
	model <- "ZIGP"
}	
if(any(model==paste("L", 1:10, sep=""))){
	true.model <- model
	model <- "ZILP"
}

models <- c("SP", "ZISP", "GP", "ZIGP", "WP", "ZIWP", "LP", "ZILP", "IP")
modelsfull <- c("single Poisson", "zero-inflated single Poisson", "gamma Poisson", "zero-inflated gamma Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson", "independant Poisson")

if((length(model) != 1) |  sum(is.na(model)) > 0 | sum(model==models)!=1){
	cat("Invalid model selection.  Please choose from ONE of the following models: ", sep="")
	cat(models, sep=", ")
	cat("\n")
	stop("Invalid model selection")
}

params.names <- switch(model, SP=c("mean"), ZISP=c("mean", "prob"), IP=c("mean", "variance"), LP=c("lmean", "lvariance"), ZILP=c("lmean", "lvariance", "prob"), GP=c("mean", "shape"), ZIGP=c("mean", "shape", "prob"), WP=c("a", "b"), ZIWP=c("a", "b", "prob"))

n.params <- switch(model, SP=1, ZISP=2, IP=2, LP=2, ZILP=3, GP=2, ZIGP=3, WP=2, ZIWP=3) 

errorpad <- quantile(0, probs=c(0.025, 0.5, 0.975, 1))
names(errorpad)[4] <- "MAX.OBS"
errorpad[] <- NA
errorpadding <- c(mean=errorpad[1], mean=errorpad[2], mean=errorpad[3], coeff.variation=errorpad[1], coeff.variation=errorpad[2], coeff.variation=errorpad[3])

if(sum(strsplit(model, split="")[[1]][1:2] != c("Z", "I"))==0){
	errorpadding <- c(errorpadding, zi=errorpad[1], zi=errorpad[2], zi=errorpad[3])
}
if((model=="GP") | (model=="WP") | (model=="ZIGP") | (model=="ZIWP")){
	errorpadding <- c(errorpadding, scale=errorpad[1], scale=errorpad[2], scale=errorpad[3], shape=errorpad[1], shape=errorpad[2], shape=errorpad[3])
}
if((model=="LP") | (model=="ZILP")){
	errorpadding <- c(errorpadding, lmean=errorpad[1], lmean=errorpad[2], lmean=errorpad[3], lvar=errorpad[1], lvar=errorpad[2], lvar=errorpad[3])
}
errorpadding <- c(errorpadding, mpsrf=NA)
if(likelihood==TRUE){
	errorpadding <- c(errorpadding, likelihood=errorpad[1], likelihood=errorpad[2], likelihood=errorpad[3], likelihood=errorpad[4])
}

if(sum(strsplit(model, split="")[[1]][1:2] != c("Z", "I"))==0){
	zero.inflation <- TRUE
}

#if(model=="SP" | model=="ZISP"){
#	n.params.inc.lambda <- n.params #  FOR PROBPOS:  + (as.numeric(likelihood) * (as.numeric(zero.inflation) * length(data)))
#}else{
#	n.params.inc.lambda <- n.params #  NOT USING LAMBDA NOW + (as.numeric(likelihood) ) # FOR PROBPOS: * (length(data) + (as.numeric(zero.inflation) * length(data))))
#}

if(model=="IP" && likelihood==TRUE){
	#  The IP model is the only one that has something monitored for lambda
	n.params.inc.lambda <- n.params + length(data)
	params.names <- c(params.names, paste("Lambda", 1:length(data)))
}else{
	n.params.inc.lambda <- n.params
}

likedata <- data
unavailablelike <- FALSE
likeadjust <- 1

if(likelihood==TRUE){
	
	if(class(data)=="list"){
		warning("The likelihood computation is not available for data provided as a list so has not been calculated")
		likelihood <- FALSE
		unavaiablelike <- TRUE
	}
	
	if(class(data)=="matrix"){
		if(ncol(data)==1){
			data <- as.integer(data)
		}else{
			if(sum(is.na(data))==0){
				likedata <- apply(data, 1, sum)
				#warning("The likelihood computation was calculated using the sum of the number of eggs calculated per animal; however this should not affect the results")
				if(model=="WP" | model=="ZIWP") warning("The likelihood computation is not available for the (ZI)WP model with repeated measures so has not been calculated")
				likeadjust <- ncol(data)
			}else{
				warning("The likelihood computation is not available for unequal repeated measures so has not been calculated")
				likelihood <- FALSE
				unavailablelike <- TRUE
			}
		}
	}
	if(likelihood) params.names <- c(params.names, "likelihood")
}	

oldtwo = oldone <- matrix(NA, ncol=n.params.inc.lambda, nrow=0)

# true.model for model testing:
strings <- run.model(model=true.model, data=data, alt.prior=alt.prior, call.jags=FALSE, monitor.lambda=likelihood, monitor.deviance=likelihood) # Only using lambda for IP model - run.model corrects for this
modelstring <- strings$model
datastring <- strings$data
new.inits <- strings$end.state
monitors <- strings$monitor

arguments$model <- modelstring
arguments$n.chains <- length(new.inits)
arguments$data <- datastring
arguments$monitor <- monitors
arguments$inits <- new.inits
arguments$plots <- TRUE
arguments$plot.type <- c('trace','density')

class(arguments) <- "list"

cat("\n")
cat("Analysing dataset using the ", model, " model...\n", sep="")
start.time <- Sys.time()

success <- try({
	#output <- autorun.jags(model=modelstring, n.chains=length(new.inits), data = datastring, monitor=monitors, silent.jags = list(silent.jags=silent.jags, killautocorr=TRUE), interactive=interactive, max.time=max.time, inits=new.inits, plots = FALSE, ...)  # OLD WAY OF CALLING SUPERSEDED BY ARGUMENTS
	output <- do.call(autorun.jags, arguments, quote=FALSE)
	}, silent=FALSE)

time.taken <- timestring(start.time, Sys.time(), units='s', show.units=FALSE)
if(class(success)=="try-error" | success[1]=="Error"){
	cat("There was an error during the simulation\n")
	errorcode <- 5
	
	if(class(success)!="try-error"){
		if(success[2] == "The simulation was aborted due to crashes"){
			cat("\nSimulation aborted due to crashes\n\n")
			errorcode <- 1
		}
	}
	
	if(raw.output==TRUE){
		return(list(mcmc=matrix(NA, ncol=n.params.inc.lambda+likelihood, nrow=1, dimnames=list(NA, params.names)), end.state=NA, samples=NA, samples.to.conv=NA, summary=NA, psrf=NA, autocorr=NA, trace=NA, density=NA, time.taken=time.taken))
	}else{
		return(c(converged=NA, error.code=errorcode, samples=NA, samples.to.conv=NA, errorpadding, time.taken=time.taken))
	}
}

samples <- output$sample
samples.to.conv <- output$burnin

psrf.target <- output$psrf$psrf.target
convergence <- output$psrf
if(any(convergence$psrf[,2]>psrf.target)){
	errorcode = 2
	converged <- FALSE
}else{
	errorcode = 0
	converged <- TRUE
}

trace <- output$trace
density <- output$density


one <- output$mcmc[[1]]
two <- output$mcmc[[2]]


if(raw.output==TRUE && likelihood==FALSE && unavailablelike==FALSE){
	if(converged==FALSE){
		cat("Returning UNCONVERGED results\n")
	}else{
		cat("Returning results\n")
	}
	cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n\n")
	return(list(mcmc=output$mcmc, end.state=output$end.state, samples=output$sample, samples.to.conv=output$burnin, summary=output$summary, psrf=output$psrf, autocorr=output$autocorr, time.taken=time.taken))
}

if(FALSE){  # NOT USING BECAUSE AUTORUN DOESN'T RETURN PARTIALLY COMPLETED CHAINS.  error and crashed don't exist - errorcode does
if(error==1) crashed <- 2

if(updates[i] >= 1000){
	if(length(one[,1]) < (1000)){
		cat("The model crashed before 1000 sampled iterations\n")
		return(c(converged=!unconverged, error.code=crashed, samples=achieved, errorpadding, time.taken=time.taken))
	}
}else{
	if(length(one[,1]) < sample){
		cat("The model crashed before ", sample, " sampled iterations\n", sep="")
		return(c(converged=!unconverged, error.code=crashed, samples=achieved, errorpadding, time.taken=time.taken))
	}
}
}  ### END DONT RUN

if(record.chains){
	chains <- output
	savingname <- paste(record.name, ".", record.setname, sep="")
	savingname <- paste(strsplit(savingname, split=" ")[[1]], collapse="_") # replace " " with "_"
	savingname <- new_unique(savingname, ".Rsave", ask=FALSE)
	save(chains, file=savingname)
}

###############  Analyse coda files

cat("Calculating results\n")

success <- try({
if(sum(strsplit(model, split="")[[1]][1:2] != c("Z", "I"))==0){
	zero.inflation <- TRUE
	model <- switch(model, ZISP="SP", ZILP="LP", ZIGP="GP", ZIWP="WP", ZILP="LP", ZIIP="IP")
}else{
	zero.inflation <- FALSE
}

shape = scale <- NA

lambda <- matrix(ncol=length(data), nrow=niter(output$mcmc)*2)

if(model=="IP"){
	mean <- c(as.matrix(one[,1]), as.matrix(two[,1]))
	variance <- c(as.matrix(one[,2]), as.matrix(two[,2]))^2
	dispersion <- variance^0.5 / mean
	if(likelihood==TRUE){
		for(i in 1:length(data)){
			lambda[,i] <- c(as.matrix(one[,i+2+as.numeric(zero.inflation)]), as.matrix(two[,i+2+as.numeric(zero.inflation)]))
			#probpos[,i] <- c(as.matrix(one[,i+length(data)+2+as.numeric(zero.inflation)]), as.matrix(two[,i+length(data)+2+as.numeric(zero.inflation)]))
		}
	}
}

if(model=="SP"){
	mean <- c(as.matrix(one[,1]), as.matrix(two[,1]))
	variance <- c(0,0)
	dispersion <- c(0,0)
	#if(likelihood==TRUE){
	#	lambda <- data
	#}
}

if(model=="LP"){
	lmean <- c(as.matrix((one[,1])), as.matrix((two[,1])))
	lvariance <- 1/c(as.matrix((one[,2])), as.matrix((two[,2])))
	normal.meanvar <- normal.params(lmean, sqrt(lvariance))
	mean <- normal.meanvar[[1]]
	variance <- normal.meanvar[[2]]^2
	dispersion <- ((exp(lvariance)) - 1)^0.5
	
	#if(exists("zilp.type")){
	#	if(zilp.type==2){
	#		lvariance <- (c(as.matrix((one[,2])), as.matrix((two[,2]))))^2
	#		normal.meanvar <- normal.params(lmean, sqrt(lvariance))
	#		mean <- normal.meanvar[[1]]
	#		variance <- normal.meanvar[[2]]^2
	#		dispersion <- ((exp(lvariance)) - 1)^0.5
	#	}
	#	if(zilp.type==3){
	#		mean <- c(as.matrix((one[,1])), as.matrix((two[,1])))
	#		lmean <- lnormal.params(mean, coeff.var=dispersion)[[1]]
	#		variance <- (normal.params(lmean, sqrt(lvariance))[[2]])^2
	#	}
	#	if(zilp.type==4){
	#		mean <- c(as.matrix((one[,1])), as.matrix((two[,1])))
	#		lvariance <- (c(as.matrix((one[,2])), as.matrix((two[,2]))))^2
	#		dispersion <- ((exp(lvariance)) - 1)^0.5
	#		lmean <- lnormal.params(mean, coeff.var=dispersion)[[1]]
	#		variance <- (normal.params(lmean, sqrt(lvariance))[[2]])^2
	#	}
	#}
	
	#if(likelihood==TRUE){
		#for(i in 1:length(data)){
			#lambda[,i] <- c(as.matrix(one[,i+2+as.numeric(zero.inflation)]), as.matrix(two[,i+2+as.numeric(zero.inflation)]))
			#probpos[,i] <- c(as.matrix(one[,i+length(data)+2+as.numeric(zero.inflation)]), as.matrix(two[,i+length(data)+2+as.numeric(zero.inflation)]))
		#}
	#}
}

if(model=="GP"){
	mean <- c(as.matrix(one[,1]), as.matrix(two[,1]))
	shape <- c(as.matrix(one[,2]), as.matrix(two[,2]))
	scale <- mean / shape
	variance <- shape*scale^2
	dispersion <- (1 / shape)^0.5
	#if(likelihood==TRUE){
	#	for(i in 1:length(data)){
	#		lambda[,i] <- c(as.matrix(one[,i+2+as.numeric(zero.inflation)]), as.matrix(two[,i+2+as.numeric(zero.inflation)])) * mean
			#probpos[,i] <- c(as.matrix(one[,i+length(data)+2+as.numeric(zero.inflation)]), as.matrix(two[,i+length(data)+2+as.numeric(zero.inflation)]))
	#	}
	#}
}

if(model=="WP"){
	#  b is different in JAGS than R, but translation done in model to make priors easier:
	# nb <- exp(-(log(b) * a)); where nb = jags 'b' and b = r 'b' (jags 'a' = r 'a')

	a <- c(as.matrix(one[,1]), as.matrix(two[,1]))
	b <- c(as.matrix(one[,2]), as.matrix(two[,2]))
	mean <- b/a * gamma(1/a)
	variance <- b^2/a * (2 * gamma(2/a) - (1 / a * gamma(1/a)^2))
	shape <- a
	scale <- b
	
	
	dispersion <- variance^0.5 / mean
	
	#if(likelihood==TRUE){
	#	for(i in 1:length(data)){
	#		lambda[,i] <- c(as.matrix(one[,i+2+as.numeric(zero.inflation)]), as.matrix(two[,i+2+as.numeric(zero.inflation)]))
			#probpos[,i] <- c(as.matrix(one[,i+length(data)+2+as.numeric(zero.inflation)]), as.matrix(two[,i+length(data)+2+as.numeric(zero.inflation)]))
	#	}
	#}
}

if(zero.inflation==TRUE & model=="SP"){
	prob <- c(as.matrix(one[,2]), as.matrix(two[,2]))
	zi <- (1 - prob) * 100
}
if(zero.inflation==TRUE & model!="SP"){
	prob <- c(as.matrix(one[,3]), as.matrix(two[,3]))
	zi <- (1- prob) * 100
}
		
if(adjust.mean==TRUE && zero.inflation==TRUE){
	mean <- mean * prob
}	
}, silent = FALSE)

if(class(success)=="try-error"){
	cat("There was an error calculating the results\n")
	if(raw.output==TRUE){
		return(list(mcmc=matrix(NA, ncol=n.params.inc.lambda+likelihood, nrow=1, dimnames=list(NA, params.names)), end.state=NA, samples=NA, samples.to.conv=NA, summary=NA, psrf=NA, autocorr=NA, trace=NA, density=NA, time.taken=time.taken))
	}else{
		return(c(converged=NA, error.code=1, samples=NA, samples.to.conv=NA, errorpadding, time.taken=time.taken))
	}
}

if(FALSE){ # TO USE quantile RATHER THAN HPDinterval
	try(zi.an <- quantile((zi), probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
	try(scale.an <- quantile(scale, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
	try(shape.an <- quantile(shape, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
	try(mean.an <- quantile(mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
	try(variance.an <- quantile(variance, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
	try(lmean.an <- quantile(lmean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
	try(lvariance.an <- quantile(lvariance, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
	try(dispersion.an <- quantile(dispersion, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
}else{ # TO USE HPDinterval
	try({hpd <- HPDinterval(coda::as.mcmc(zi))
		zi.an <- c(hpd[1], median(zi), hpd[2])
		names(zi.an) <- c("l.95", "median", "u.95")}, silent=TRUE)
	try({hpd <- HPDinterval(coda::as.mcmc(scale))
		scale.an <- c(hpd[1], median(scale), hpd[2])
		names(scale.an) <- c("l.95", "median", "u.95")}, silent=TRUE)
	try({hpd <- HPDinterval(coda::as.mcmc(shape))
		shape.an <- c(hpd[1], median(shape), hpd[2])
		names(shape.an) <- c("l.95", "median", "u.95")}, silent=TRUE)
	try({hpd <- HPDinterval(coda::as.mcmc(mean))
		mean.an <- c(hpd[1], median(mean), hpd[2])
		names(mean.an) <- c("l.95", "median", "u.95")}, silent=TRUE)
	try({hpd <- HPDinterval(coda::as.mcmc(variance))
		variance.an <- c(hpd[1], median(variance), hpd[2])
		names(variance.an) <- c("l.95", "median", "u.95")}, silent=TRUE)
	try({hpd <- HPDinterval(coda::as.mcmc(lmean))
		lmean.an <- c(hpd[1], median(lmean), hpd[2])
		names(lmean.an) <- c("l.95", "median", "u.95")}, silent=TRUE)
	try({hpd <- HPDinterval(coda::as.mcmc(lvariance))
		lvariance.an <- c(hpd[1], median(lvariance), hpd[2])
		names(lvariance.an) <- c("l.95", "median", "u.95")}, silent=TRUE)
	try({hpd <- HPDinterval(coda::as.mcmc(dispersion))
		dispersion.an <- c(hpd[1], median(dispersion), hpd[2])
		names(dispersion.an) <- c("l.95", "median", "u.95")}, silent=TRUE)
}

autocorrdis <- NA
try({
	autocorrdis <- coda::autocorr.diag(coda::as.mcmc(dispersion))[2,1]
	names(autocorrdis) <- NULL
	}, silent=TRUE)

results <- c(converged=as.integer(converged), error.code=as.integer(errorcode), samples=as.integer(samples), samples.to.conv=as.integer(samples.to.conv), mean=mean.an[1], mean=mean.an[2], mean=mean.an[3], coeff.variation=dispersion.an[1], coeff.variation=dispersion.an[2], coeff.variation=dispersion.an[3])

if(zero.inflation==TRUE){
	results <- c(results, zi=zi.an[1], zi=zi.an[2], zi=zi.an[3])
}

if((model=="GP") | (model=="WP")){
	results <- c(results, scale=scale.an[1], scale=scale.an[2], scale=scale.an[3], shape=shape.an[1], shape=shape.an[2], shape=shape.an[3])
}
if(model=="LP"){
	results <- c(results, lmean=lmean.an[1], lmean=lmean.an[2], lmean=lmean.an[3], lvar=lvariance.an[1], lvar=lvariance.an[2], lvar=lvariance.an[3])
}

if(n.params==1) mpsrf <- convergence$psrf[1,1] else mpsrf <- convergence$mpsrf

if(class(mpsrf)=="character") mpsrf <- NA

results <- c(results, mpsrf=mpsrf)

if(likelihood==TRUE){
	# Old likelihoods calculated manually:
	if(TRUE){
		
	suppressWarnings(success <- try({
	if(zero.inflation==FALSE){
		l.zi <- NA
	}else{
		l.zi <- zi
	}

	if(model=="WP" | model=="GP"){
		# likeadjust has to be 1 for (ZI)WP model
		likeli <- likelihood(model=model, data=likedata, shape=shape, scale=scale*likeadjust, zi=l.zi, silent=paste("noziwarn", silent.jags, sep=""), log=TRUE, raw.output=TRUE)  # Leave iterations as default (1000)
	}
	if(model=="LP"){
		newmeanvar <- normal.params(lmean, sqrt(lvariance))
		likemean <- newmeanvar[[1]] * likeadjust
		likecoeff <- newmeanvar[[3]]
		likevar <- (likemean/likecoeff)^2
		likeli <- likelihood(model=model, data=likedata, mean=likemean, variance=likevar, zi=l.zi, silent=paste("noziwarn", silent.jags, sep=""), log=TRUE, raw.output=TRUE)  # Leave iterations as default (1000)
		
	}
	if(model=="SP"){
		likeli <- likelihood(model=model, data=likedata, mean=mean*likeadjust, zi=l.zi, silent=paste("noziwarn", silent.jags, sep=""), log=TRUE, raw.output=TRUE)  # Leave iterations as default
	}
	if(model=="IP"){
		likeli <- likelihood(model=model, data=likedata, mean=lambda*likeadjust, zi=l.zi, silent=paste("noziwarn", silent.jags, sep=""), log=TRUE, raw.output=TRUE)
	}}, silent=FALSE))
	
	if(class(success)=="try-error"){
		cat("An error occured while computing the likelihood\n")
		likeli <- list(NA, NA)
	}else{
		if(silent.jags) cat("Finished calculating the likelihood\n")
	}
	
	}else{
		# New likelihoods based on deviance:
		likeli <- as.numeric(combine.mcmc(output$mcmc)[,"deviance"])/-2
		
	}
	
	#print(sum(is.na(likeli[[1]])))
	#assign('like', likeli, pos=.GlobalEnv)
	if(any(is.na(likeli[[1]]))){
		likeli.an <- c(l.95=NA, median=NA, u.95=NA, MAX.OBS=NA)
	}else{
		hpd <- HPDinterval(coda::as.mcmc(likeli[[1]]))
		likeli.an <- c(l.95=hpd[1], median=median(likeli[[1]]), u.95=hpd[2], MAX.OBS=max(likeli[[1]]))
	}
	results <- c(results, likelihood=likeli.an[1], likelihood=likeli.an[2], likelihood=likeli.an[3], likelihood=likeli.an[4])
}

if(unavailablelike==TRUE){
	likeli.an <- quantile(NA, probs=c(0.025, 0.5, 0.975, 1), na.rm=TRUE)
	names(likeli.an)[4] <- "MAX.OBS"
	results <- c(results, likelihood=likeli.an[1], likelihood=likeli.an[2], likelihood=likeli.an[3], likelihood=likeli.an[4])
}

if(raw.output==TRUE & (likelihood==TRUE | unavailablelike==TRUE)){
	
	# No longer need to append likelihoods to each chain - deviance is there already:
	#if(FALSE){
	success <- try({
	
	one.backup <- one
	two.backup <- two
	
	newone <- matrix(NA, nrow=nrow(one), ncol=ncol(one)+1)
	newtwo <- matrix(NA, nrow=nrow(two), ncol=ncol(two)+1)
	
	newone[,1:ncol(one)] <- one
	newtwo[,1:ncol(two)] <- two
	
	dimnames(newone) <- list(NULL, c(dimnames(one)[[2]], "likelihood"))
	dimnames(newtwo) <- list(NULL, c(dimnames(two)[[2]], "likelihood"))	
	
	one <- newone
	two <- newtwo
	
	sequence <- numeric(length=(length(one[,1]) + length(two[,1])))
	
	sequence[] <- NA
	
	if(unavailablelike==FALSE){
		sequence[likeli[[2]]] <- likeli[[1]]
	}	
	
	one[,length(one[1,])] <- sequence[1:length(one[,1])]
	two[,length(two[1,])] <- sequence[(length(one[,1])+1):(length(one[,1])+length(two[,1]))]
	
	}, silent=FALSE)
	
	if(class(success)=="try-error"){
		cat("An error occured while combining the chains with the calculated likelihoods.  The likelihoods will not be returned\n")
		one <- one.backup
		two <- two.backup
	}
	#}
	
	mcmc <- coda::mcmc.list(coda::as.mcmc(one), coda::as.mcmc(two))

#	mcmc[[1]][,"deviance"] <- mcmc[[1]][,"deviance"]/-2
#	mcmc[[2]][,"deviance"] <- mcmc[[2]][,"deviance"]/-2
	
#	varnames(mcmc)[nvar(mcmc)] <- "likelihood"
	
	if(converged==FALSE){
		cat("Returning UNCONVERGED results\n")
	}else{
		cat("Returning results\n")
	}
	cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n\n")
	
	return(list(mcmc=mcmc, end.state=output$end.state, samples=output$sample, samples.to.conv=output$burnin, summary=output$summary, psrf=output$psrf, autocorr=output$autocorr, trace=trace, density=density, time.taken=time.taken))
}

largeod <- assess.variance(model=list(model=model, alt.prior=alt.prior, l.95 = dispersion.an[1], u.95 = dispersion.an[3]))

if((alt.prior==TRUE | is.character(alt.prior)) && largeod==TRUE){
	cat("*WARNING*  The 95% confidence interval for the variance is very large, indicating a lack of information about this parameter in the data.  The model should be re-run with the standard prior distribution and the two sets of results compared to ensure that the results for all parameters are reliable\n")
}
if(alt.prior==FALSE && largeod==TRUE){
	cat("*WARNING*  The 95% confidence interval for the variance is very large, indicating a lack of information about this parameter in the data.  The model should be re-run with an alternative prior distribution and the two sets of results compared to ensure that the results for all parameters are reliable\n")
}


cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n")
cat("Finished running the model\n\n", sep="")
# results <- c(results, autocorr.disp = autocorrdis, time.taken=time.taken)  # Removing autocorr.disp for now as not particularly useful

results <- c(results, time.taken=time.taken)
names <- names(results)
results <- as.numeric(results)
names(results) <- names
return(results)

}	

fec.analysis <- bayescount.single
FEC.analysis <- bayescount.single
count.analysis <- bayescount.single

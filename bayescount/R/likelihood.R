likelihood <- function(model=stop("Please specify a distribution"), data=stop("Data must be specified"), mean=NA, variance=NA, zi=NA, shape=NA, scale=NA, iterations=min(1000, (length(as.matrix(mean)[,1])*length(shape))), log=TRUE, silent=FALSE, raw.output=FALSE){
	
	model <- toupper(model)
	

	## Cheat to stop function warning of non-zero inflated model when called from bayescount.single() and multiple deprecation warnings
	zi.warn <- TRUE
	if(class(silent)=='logical')
		warning('The likelihood function is deprecated and will be removed from version 1.0 of bayescount (expected to be released in mid 2015)')
	if(silent=='TRUE')
		silent <- TRUE
	if(silent=='FALSE')
		silent <- FALSE
	
	if(silent=="noziwarnTRUE"){
		zi.warn <- FALSE
		silent <- TRUE
	}
	if(silent=="noziwarnFALSE"){
		zi.warn <- FALSE
		silent <- FALSE
	}


	## cheat to allow IP model to use SP model stuff:
	
	if(model=="SP" | model=="ZISP" | model=="P" | model=="ZIP"){
		if(!is.matrix(mean)){
			oldmean <- mean
			mean <- matrix(NA, ncol=length(na.omit(data)), nrow=length(oldmean))
			mean[] <- oldmean
		}
	}
	
	model <- switch(model, SP="P", IP="P", ZISP="ZIP", model)
	models <- c("P", "ZIP", "G", "ZIG", "L", "ZIL", "W", "ZIW", "GP", "ZIGP", "LP", "ZILP", "WP", "ZIWP")
	
	if((length(model) != 1) |  sum(is.na(model)) > 0 | sum(model==models)!=1){
		if(silent==FALSE){
			cat("Invalid model selection.  Please choose from ONE of the following distributions: ", sep="")
			cat(models, sep=", ")
			cat("\n")
		}
		stop("Invalid model selection")
	}
	
	suppressWarnings(if((sum(na.omit(strsplit(model, split=""))[[1]][] == c("Z", "I")))==2){
		if(is.na(zi)){
			stop("A value for zero-inflation (%) must be supplied when using a zero-inflated model")
		}
		zero.inflation <- TRUE
		prob <- 1 - (zi/100)
		model <- switch(model, ZIP="P", ZIG="G", ZIW="W", ZIL="L", ZIGP="GP", ZILP="LP", ZIWP="WP")
	}else{
		if(!is.na(zi)){
			if(zi.warn==TRUE){
				cat("A value for zero-inflation was provided for a non zero-inflated model.  The appropriate zero-inflated model will be used instead of the model specified\n")
			}
			zero.inflation <- TRUE
			prob <- 1 - (zi/100)
		}else{
			zero.inflation <- FALSE
		}
	}
	)
	
	
	chain.length <- length(as.matrix(mean)[,1]) * length(shape)
	
	if(chain.length < iterations){
		stop("The number of iterations to calculate the likelihood for must not be more than the length of the chain provided")
	}
	
	s.info <- Sys.info()
	p.info <- .Platform
	
	if(as.numeric(paste(R.version$major, R.version$minor, sep="")) < 26){
    	stop("Please update your version of R to the latest available (> 2.6.0 required)")
	}
	
	os <- p.info$OS.type
	username <- as.character(s.info["user"])
	rversion <- R.version$version
	gui <- p.info$GUI
	p.type <- p.info$pkgType
	
	guitest <- c("R.GUI"=gui, "R.package.type"=p.type, "R.version"=list(rversion))
	
	#if(guitest$os=='windows' | (guitest$R.GUI == "AQUA" & guitest$R.package.type == "mac.binary")) macgui <- TRUE else macgui <- FALSE
	
	if(guitest$R.GUI == "AQUA" & guitest$R.package.type == "mac.binary") macgui <- TRUE else macgui <- FALSE
	
	###  The GUI version of R on mac currently doesn't support destructive backspace or carriage return, so the following is an alternative way of watching progress for Mac GUI
	##  later discovered that Windows doesn't handle them properly either - not sure about this now
	
	sequence <- pmin(round((1:iterations) * (chain.length / iterations), digits=0), chain.length)
	iters <- length(sequence)	
	progs <- character(length=iters)
	progs[] <- ""
	
	if(sum(is.na(data)) > 0){
		if(silent==FALSE) cat("WARNING:  Missing data was removed before calculating the likelihood\n")
		data <- na.omit(data)
	}
	
	results <- matrix(NA, ncol=length(data), nrow=iters)
	
	if(silent==FALSE){
		cat("Calculating likelihood for the ", if(zero.inflation) "ZI", model, " model\n", sep="")
		if(macgui==TRUE){
			on.exit(cat("Finished calculating the likelihood\n"))
		}else{
			on.exit(cat("\b\rFinished calculating the likelihood                    \n"))
		}
	}
	
	limits.data <- data
	limits.data[limits.data==0] <- 1
	
	integrate.min <- qgamma(0.001, shape=limits.data, rate=1)
	integrate.max <- qgamma(0.999, shape=limits.data, rate=1)
	
	integrate.min[data==0] <- 0
	
	if(macgui==TRUE){
		
		progress.function <- function(iteration){
			cat(progs[iteration])
		}
		
		if(silent==FALSE && (model=="LP" | model=="GP" | model=="WP")){
			cat("0%---------------------50%--------------------100%\n")
			if(iterations >=50){
				progs[pmin(round((1:50)*(iters/50), digits=0), iters)] <- "*"
			}
			if(iterations < 50 && iterations >= 25){	
				progs[pmin(round((1:25)*(iters/25), digits=0), iters)] <- "**"
			}
			if(iterations < 25 && iterations >= 10){	
				progs[pmin(round((1:10)*(iters/10), digits=0), iters)] <- "*****"
			}
			if(iterations < 10 && iterations >= 5){	
				progs[pmin(round((1:5)*(iters/5), digits=0), iters)] <- "**********"
			}
			if(iterations == 4){	
				progs[c(1,3)] <- "************"
				progs[c(2,4)] <- "*************"
			}
			if(iterations == 3){	
				progs[c(1,3)] <- "*****************"
				progs[2] <- "****************"
			}
			if(iterations == 2){	
				progs[] <- "*************************"
			}
			if(iterations == 1){	
				progs <- "**************************************************"
			}
			progs[iters] <- paste(progs[iters], "\n", sep="")
		}
		
	}else{
	
		start.time <- Sys.time()
	
		progress.function <- function(iteration){
			
			post.time <- Sys.time()
			total.time <- timestring(start.time, post.time, show.units=TRUE)
			percent.complete <- iteration / iters
			time.remaining <- timestring((as.integer(difftime(post.time, start.time, units="secs")) / percent.complete) - as.integer(difftime(post.time, start.time, units="secs")), show.units=TRUE)
			if(iteration!=1){
				cat("\b\r")
			}
			cat(round(percent.complete*100), "% complete, approximately ", time.remaining, " remaining      ", sep="")
		
		}
	}
	
	if(silent==TRUE){
		progress.function <- function(iteration){
		}
	}
	
	
	if(model=="P"){
		if(sum(is.na(mean)) > 0){
			stop("A value for mean must be supplied when using the Poisson model")
		}
		if(sum(!is.na(variance)) > 0 && silent==FALSE){
			cat("The value for variance provided will be ignored using the Poisson model\n")
		}
		if(sum(!is.na(shape)) > 0 && silent==FALSE){
			cat("The value for the shape parameter provided will be ignored using the Poisson model\n")
		}
		if(sum(!is.na(scale)) > 0 && silent==FALSE){
			cat("The value for the scale parameter provided will be ignored using the Poisson model\n")
		}
		
		if(zero.inflation==TRUE){
						
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					results[iteration,datapoint] <- log((((data[datapoint]==0) * (1 - prob[sequence[iteration]])) + (prob[sequence[iteration]] * exp(-mean[sequence[iteration], datapoint])*mean[sequence[iteration], datapoint]^data[datapoint] / factorial(data[datapoint]))))
				
				}
			}
			
		}else{
			
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					results[iteration,datapoint] <- log(dpois(data[datapoint], mean[sequence[iteration], datapoint]))
				
				}
			}
		}
		
		likelihoods <- apply(results, 1, sum)
		
	}
	
	if(model=="L"){
		if(sum(is.na(mean))>0){
			stop("A value for mean must be supplied when using the lognormal model")
		}
		if(sum(is.na(variance))>0){
			stop("A value for variance must be supplied when using the lognormal model")
		}
		if(sum(!is.na(shape))>0 && silent==FALSE){
			cat("The value for the shape parameter provided will be ignored using the lognormal model\n")
		}
		if(sum(!is.na(scale))>0 && silent==FALSE){
			cat("The value for the scale parameter provided will be ignored using the lognormal model\n")
		}
		
		
		newlnorms <- lnormal.params(mean, sqrt(variance))
		lmean <- newlnorms[[1]]
		lsd <- newlnorms[[2]]
				
		if(zero.inflation==TRUE){
		
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					results[iteration,datapoint] <- log(((data[datapoint]==0) * (1 - prob[sequence[iteration]])) + (prob[sequence[iteration]] * dlnorm(data[datapoint], lmean[sequence[iteration]], lsd[sequence[iteration]])))
				
				}
			}
		
		}else{
			
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
							
					results[iteration,datapoint] <- log(dlnorm(data[datapoint], lmean[sequence[iteration]], lsd[sequence[iteration]]))
				
				}
			}
		
		}
		
		likelihoods <- apply(results, 1, sum)		
	}
	
	if(model=="LP"){
		if(any(is.na(mean))){
			stop("A value for mean must be supplied when using the lognormal Poisson model")
		}
		if(any(is.na(variance))){
			stop("A value for variance must be supplied when using the lognormal Poisson model")
		}
		if(any(!is.na(shape)) && silent==FALSE){
			cat("The value for the shape parameter provided will be ignored using the lognormal Poisson model\n")
		}
		if(any(!is.na(scale)) && silent==FALSE){
			cat("The value for the scale parameter provided will be ignored using the lognormal Poisson model\n")
		}
		
		newlnorms <- lnormal.params(mean, sqrt(variance))
		lmean <- newlnorms[[1]]
		lsd <- newlnorms[[2]]
		
		suppressWarnings({
		
		if(zero.inflation==TRUE){
		
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					f <- function(lambda) dpois(data[datapoint], lambda) * dlnorm(lambda, lmean[sequence[iteration]], lsd[sequence[iteration]])
									
					results[iteration,datapoint] <- try(log(((data[datapoint]==0) * (1 - prob[sequence[iteration]])) + (prob[sequence[iteration]] * integrate(f, integrate.min[datapoint], integrate.max[datapoint], stop.on.error=FALSE)$value)), silent=TRUE)
				
				}
			}
		
		}else{
			
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					f <- function(lambda) dpois(data[datapoint], lambda) * dlnorm(lambda, lmean[sequence[iteration]], lsd[sequence[iteration]])
			
					results[iteration,datapoint] <- try(log(integrate(f,integrate.min[datapoint], integrate.max[datapoint], stop.on.error=FALSE)$value), silent=TRUE)
					
				
				}
			}
		
		}
		results <- matrix(data=as.numeric(results), nrow=length(results[,1]), ncol=length(results[1,]))})
		
		likelihoods <- apply(results, 1, sum)
}
	
	if(model=="G"){
		if(any(is.na(shape))){
			stop("A value for the shape parameter must be supplied when using the gamma model")
		}
		if(any(is.na(scale))){
			stop("A value for the scale parameter must be supplied when using the gamma model")
		}
		if(any(!is.na(mean)) && silent==FALSE){
			cat("The value for mean provided will be ignored using the gamma model\n")
		}
		if(any(!is.na(variance)) && silent==FALSE){
			cat("The value for variance provided will be ignored using the gamma model\n")
		}
		
		if(zero.inflation==TRUE){
		
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					results[iteration,datapoint] <- log(((data[datapoint]==0) * (1 - prob[sequence[iteration]])) + (prob[sequence[iteration]] * dgamma(data[datapoint], shape=shape[sequence[iteration]], scale=scale[sequence[iteration]])))
				
				}
			}
		
		}else{
			
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
							
					results[iteration,datapoint] <- log(dgamma(data[datapoint], shape=shape[sequence[iteration]], scale=scale[sequence[iteration]]))
				
				}
			}
		
		}
		
		likelihoods <- apply(results, 1, sum)
		
}
	
	if(model=="GP"){
		if(any(is.na(shape))){
			stop("A value for the shape parameter must be supplied when using the gamma Poisson model")
		}
		if(any(is.na(scale))){
			stop("A value for the scale parameter must be supplied when using the gamma Poisson model")
		}
		if(any(!is.na(mean)) && silent==FALSE){
			cat("The value for mean provided will be ignored using the gamma Poisson model\n")
		}
		if(any(!is.na(variance)) && silent==FALSE){
			cat("The value for variance provided will be ignored using the gamma Poisson model\n")
		}
		
		suppressWarnings({
		
		if(zero.inflation==TRUE){
		
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					f <- function(lambda) dpois(data[datapoint], lambda) * dgamma(lambda, shape=shape[sequence[iteration]], scale=scale[sequence[iteration]])
			
					results[iteration,datapoint] <- try(log(((data[datapoint]==0) * (1 - prob[sequence[iteration]])) + (prob[sequence[iteration]] * integrate(f, integrate.min[datapoint], integrate.max[datapoint], stop.on.error=FALSE)$value)), silent=TRUE)
				
				}
			}
		
		}else{
			
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					f <- function(lambda) dpois(data[datapoint], lambda) * dgamma(lambda, shape=shape[sequence[iteration]], scale=scale[sequence[iteration]])
			
					results[iteration,datapoint] <- try(log(integrate(f, integrate.min[datapoint], integrate.max[datapoint], stop.on.error=FALSE)$value), silent=TRUE)
				
				}
			}
		
		}
		results <- matrix(data=as.numeric(results), nrow=length(results[,1]), ncol=length(results[1,]))})
		
		likelihoods <- apply(results, 1, sum)
		
}
	
	if(model=="W"){
		if(any(is.na(shape))){
			stop("A value for the shape parameter must be supplied when using the Weibull model")
		}
		if(any(is.na(scale))){
			stop("A value for the scale parameter must be supplied when using the Weibull model")
		}
		if(any(!is.na(mean)) && silent==FALSE){
			cat("The value for mean provided will be ignored using the Weibull model\n")
		}
		if(any(!is.na(variance)) && silent==FALSE){
			cat("The value for variance provided will be ignored using the Weibull model\n")
		}
		
		if(zero.inflation==TRUE){
		
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					results[iteration,datapoint] <- log(((data[datapoint]==0) * (1 - prob[sequence[iteration]])) + (prob[sequence[iteration]] * dweibull(data[datapoint], shape=shape[sequence[iteration]], scale=scale[sequence[iteration]])))
				
				}
			}
		
		}else{
			
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
							
					results[iteration,datapoint] <- log(dweibull(data[datapoint], shape=shape[sequence[iteration]], scale=scale[sequence[iteration]]))
				
				}
			}
		
		}
		
		likelihoods <- apply(results, 1, sum)
		
}
	
	if(model=="WP"){
		if(any(is.na(shape))){
			stop("A value for the shape parameter must be supplied when using the Weibull Poisson model")
		}
		if(any(is.na(scale))){
			stop("A value for the scale parameter must be supplied when using the Weibull Poisson model")
		}
		if(any(!is.na(mean)) && silent==FALSE){
			cat("The value for mean provided will be ignored using the Weibull Poisson model\n")
		}
		if(any(!is.na(variance)) && silent==FALSE){
			cat("The value for variance provided will be ignored using the Weibull Poisson model\n")
		}
		
		suppressWarnings({
		
		if(zero.inflation==TRUE){
		
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					f <- function(lambda) dpois(data[datapoint], lambda) * dweibull(lambda, shape=shape[sequence[iteration]], scale=scale[sequence[iteration]])
			
					results[iteration,datapoint] <- try(log(((data[datapoint]==0) * (1 - prob[sequence[iteration]])) + (prob[sequence[iteration]] * integrate(f, integrate.min[datapoint], integrate.max[datapoint], stop.on.error=FALSE)$value)), silent=TRUE)
				
				}
			}
		
		}else{
			
			for(iteration in 1:iters){
				progress.function(iteration)
				for(datapoint in 1:length(data)){
					
					f <- function(lambda) dpois(data[datapoint], lambda) * dweibull(lambda, shape=shape[sequence[iteration]], scale=scale[sequence[iteration]])
			
					results[iteration,datapoint] <- try(log(integrate(f, integrate.min[datapoint], integrate.max[datapoint], stop.on.error=FALSE)$value), silent=TRUE)
				
				}
			}
		
		}
		results <- matrix(data=as.numeric(results), nrow=length(results[,1]), ncol=length(results[1,]))})
#		assign('gresults', results, pos=.GlobalEnv)
		likelihoods <- apply(results, 1, sum)
		
}
	
	
	if(log==FALSE){
		likelihoods <- exp(likelihoods)
	}
	#if(silent==FALSE) cat("\n")
	if(any(is.na(results)) && silent==FALSE && any(c(model=="WP", model=="GP", model=="LP"))) cat("Error:  The likelihood could not be calculated at one or more iterations because an infinitely small likelihood was integrated for one or more datapoints at these iterations\n")
	
	if(raw.output==TRUE){
		return(list(likelihood=likelihoods, iteration=sequence))
	}else{
		if(length(likelihoods)==1){
			return(likelihoods)
		}else{
			if(any(is.na(results))){
				likeli.an <- c(l.95=NA, median=NA, u.95=NA, MAX=NA)
			}else{
				hpd <- coda::HPDinterval(coda::as.mcmc(likelihoods))
				likeli.an <- c(l.95=hpd[1], median=median(likelihoods), u.95=hpd[2], MAX=max(likelihoods))
			}
			return(likeli.an)
		}
	}

}
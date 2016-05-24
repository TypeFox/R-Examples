fecrt.analysis <- function(name = NA, pre.data = NA, post.data = NA, data = list(pre=pre.data, post=post.data), animal.names = FALSE, efficacy=95, confidence=0.95, control.animals=FALSE, divide.data = 1, record.chains = FALSE, write.file = FALSE, bootstrap.iters=10000, plot.graph = TRUE, skip.mcmc = FALSE, ...){

runname <- name

# Individual analysis is removed for now, unless I add a distribution of efficacy for the paired model later.
# individual.analysis=paired.model, 
#if(!paired.model & individual.analysis) stop("Individual analysis is only available using the paired model")

lci <- 0+((1-(confidence))/2)
uci <- 1-((1-(confidence))/2)

div <- divide.data 

passthrough <- list(...)
# Override some defaults:
if(is.null(passthrough$max.time))
	passthrough$max.time <- "1hr"
if(is.null(passthrough$interactive))
	passthrough$interactive <- FALSE

stopifnot(length(confidence)==1 && confidence>0 && confidence<1)
passthrough$confidence <- confidence

if(any(names(passthrough)=='paired.model'))
	paired.model <- passthrough$paired.model
else
	paired.model <- FALSE

#arguments <- formals(autoextend.jags)

#for(i in 1:length(passthrough)){
#	if(is.null(arguments[[names(passthrough)[i]]])){
		#warning(paste("'", names(passthrough)[i], "' is not an argument to autoextend.jags and was ignored"))
#	}else{
#		arguments[[names(passthrough)[i]]] <- passthrough[[i]]
#	}
#}

arguments <- passthrough[names(passthrough)%in%names(formals(autoextend.jags))]

testwritable <- new_unique("test")
if(testwritable=="Directory not writable"){
	cat("\nThe working directory is not writable.  Please change the working directory\n\n")
	stop("Directory not writable")
}

if(class(data)=="function") stop("The class of the object supplied for data is a function")

if(class(data)!="character" & class(data)!="matrix"){
		
	if(class(data)!="list") stop("Data supplied in an incorrect format")
	if(length(data)!=2) stop("Data supplied in an incorrect format - the list must be of length two")
	
	pre.data <- data[[1]]
	post.data <- data[[2]]
	
	if(class(pre.data)=="integer") pre.data <- as.numeric(pre.data)
	if(class(post.data)=="integer") post.data <- as.numeric(post.data)
	
	if(class(pre.data)!=class(post.data)) stop("The pre and post treatment data must be provided in the same format")
	
	if(class(pre.data)=="array" & paired.model==FALSE) stop("The paired model is required when using repeat samples within animal.  Either specify the data as a matrix, or set paired.model=TRUE")

	if(class(pre.data)=="matrix"){
		dims <- dim(pre.data)
		pre.data <- array(t(pre.data), dim=c(1,dims[2],dims[1]))
		dims <- dim(post.data)
		post.data <- array(t(post.data), dim=c(1,dims[2],dims[1]))
	}
	
	if(class(pre.data)=="numeric" | class(pre.data)=="integer"){
		pre.data <- array(pre.data, dim=c(1,1,length(pre.data)))
		post.data <- array(post.data, dim=c(1,1,length(post.data)))
	}
	
	if(dim(pre.data)[3] != dim(post.data)[3]) stop("Unequal numbers of animals pre- and post-treatment")
	
	data <- list(pre.counts=pre.data, post.counts=post.data)
	
}



cat("\n--- FECRT: Analyse faecal egg cout reduction test data using a Bayesian distributional simulation model ---\n\n")
cat("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*\n\n")

if(is.na(runname)){
	runname <- ask(prompt = "Please enter a name for this analysis:  ", type="character")
}
name <- runname

datain <- data
datana <- data
datana <- as.vector(na.omit(datana))
length(datana) <- 1

dataok=FALSE

if(class(data)=="matrix" | class(data)=="list"){
	dataok <- TRUE
}else{
	if(is.na(datana)==FALSE){
	exists <- try(file.exists(datana), silent=TRUE)
	if((class(exists)=="try-error")==FALSE){
		if(exists==TRUE){
			suppressWarnings(data <- try(read.csv(datain, header=FALSE), silent=TRUE))
		}
	}
	suppressWarnings(valid.data <- try((length(as.matrix(data[,1])) > 1), silent=TRUE))
	if((class(valid.data)=="try-error")==TRUE){
		cat("ERROR:  The path you have entered does not appear to be valid\n")
	}else{
		if(valid.data==FALSE){
			cat("ERROR:  Invalid path / data\n") 
		}else{
			dataok=TRUE
		}
	}
	cat("\n")
	}
}
while(dataok==FALSE){
	datain <- ask(prompt = "Please enter the path to a (comma delimited) CSV file containing the data (type 'exit' to quit):  ", type="character")
	if((datain=="exit")==TRUE){
		stop("User exited the program")
	}
	exists <- try(file.exists(datain), silent=TRUE)
	if((class(exists)=="try-error")==FALSE){
		if(exists==TRUE){
		data <- try(as.matrix(read.csv(datain, header=FALSE), silent=TRUE))
		if((class(data)=="try-error")==FALSE){
			valid.data <- try(length(data[,1]) > 1, silent=TRUE)
			if((class(valid.data)=="try-error")==TRUE){
				cat("ERROR:  The path you have entered does not appear to be valid\n")
			}else{
				if(valid.data==FALSE){
					cat("ERROR:  The data you have entered is of length less than 2\n")
				}else{
					dataok=TRUE
				}
			}
		}else{
			cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
		}
		}else{
			cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
		}
	}else{
		cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
	}
	cat("\n")

}

if(class(data)!="matrix" & identical(animal.names,TRUE)){
	warning("'animal.names' cannot be TRUE if the data is not provided as a matrix.  'animal.names' will be set to FALSE")
	animal.names <- FALSE
}

if(class(data)=="matrix"){
	
	if(ncol(data) < (2+identical(animal.names,TRUE)) | ncol(data) > (3+identical(animal.names,TRUE))) stop("If provided as a matrix, the data must have either 2, 3 or 4 columns (if animal.names==TRUE)")
	
	if(identical(animal.names,TRUE)){
		animal.names <- data[,1]
		data <- data[,2:ncol(data)]
	}
	if(ncol(data)==3) control.animals <- data[,3]

	pre.data <- array(data[,1], dim=c(1,1,length(data[,1])))
	post.data <- array(data[,2], dim=c(1,1,length(data[,2])))

	data <- list(pre.counts=pre.data, post.counts=post.data)
}

if(class(data$pre.counts) == "NULL" | class(data$post.counts) == "NULL") stop("An error occured while transforming the data to an appropriate format")

if(!identical(animal.names, FALSE)){
	if(length(animal.names) != length(data$pre.counts[1,1,])) stop("The length of the character vector animal.names does not match the length of the data provided")
	names <- animal.names
}else{
	names <- paste("Animal ", 1:length(data$pre.counts[1,1,]), sep="")
}

N <- length(data$pre.counts[1,1,])

if(N > 100 & !skip.mcmc){
	if(arguments$interactive){
		if(!speed){
			cat("There are a large number of animals in the data which will take a long time to analyse using MCMC, and may cause memory issues with individual analysis.  ")
			returned <- ask("Please choose from the following options:\n1 Continue with analysis (this may cause a crash)\n2 Set individual.analysis to FALSE\n3 Skip the MCMC analysis\n:", "numeric", c(1,3))
			cat("\n")
			if(returned==2){
				individual.analysis <- FALSE
				speed <- TRUE
			}
			if(returned==3) skip.mcmc <- TRUE
		}else{
			cat("There are a large number of animals in the data which will take a long time to analyse using MCMC.  ")
			returned <- ask("Please choose from the following options:\n1 Continue with analysis (this may take a while)\n2 Skip the MCMC analysis\n:", "numeric", c(1,2))
			cat("\n")
			if(returned==2) skip.mcmc <- TRUE
		}
	}else{
		if(!speed){
			individual.analysis <- FALSE
			speed <- TRUE
			warning("There are a large number of animals in the data which will take a long time to analyse using MCMC, and may cause memory issues with individual analysis.  individual analysis was therefore not performed.")
		}
	}
}

pre <- data$pre.counts / div
post <- data$post.counts / div

if(identical(control.animals, FALSE)){
	control.animals <- replicate(N, FALSE)
}

control.animals <- as.integer(control.animals)

if(length(control.animals)!=N) stop("The length of the control/treatment vector does not match the number of animals")

if(length(pre)!=length(post)) warning("Pre and post treatment data are of different lengths")

if(sum(control.animals)==0) cat("Assessing the faecal egg count reduction for ", N, " animals.  This will take some time...\n", sep="") else cat("Assessing the faecal egg count reduction for ", N, " animals (including ", sum(control.animals==1), " control animals).  This will take some time...\n", sep="")



# Individual analysis has been removed for now; may add it back in later if I use varying efficacy
#if(individual.analysis) monitor <- c(monitor, "ind.delta.mean") #individual.analysis can only happen for paired model
#if(zero.inflation & individual.analysis) monitor <- c(monitor, "probpos")

if(!skip.mcmc){

	# Get model from fecrt.model:
	
	prean=presamp=precont <- pre
	for(a in 1:N){
		prean[,,a] <- a
		precont[,,a] <- control.animals[a]+1
	}
	for(s in 1:dim(presamp)[2]){
		presamp[,s,] <- s
	}
	postan=postsamp=postcont <- post
	for(a in 1:N){
		postan[,,a] <- a
		postcont[,,a] <- control.animals[a]+1
	}
	for(s in 1:dim(postsamp)[2]){
		postsamp[,s,] <- s
	}
	modeldata <- data.frame(Count=c(as.numeric(pre), as.numeric(post)), Sample=c(as.numeric(presamp), as.numeric(postsamp)), Subject=c(as.numeric(prean), as.numeric(postan)), Control=c(as.numeric(precont), as.numeric(postcont)), Time=c(rep(1, length(as.numeric(pre))), rep(2, length(as.numeric(post)))))
	
	arglist <- passthrough
	arglist$data <- modeldata
	
	rjo <- do.call('fecrt.model', arglist, quote=FALSE)
	
	arguments$runjags.object <- rjo
	class(arguments) <- 'list'
	results <- do.call('autoextend.jags', arguments, quote=FALSE)
	
	#results <- autorun.jags(data=datastring, model=model, monitor=monitor, n.chains=2, inits=c(inits1, inits2), silent.jags = list(silent.jags=silent.jags, killautocorr=TRUE), plots = FALSE, thin.sample = TRUE, interactive=interactive, max.time=max.time, ...)

if(results[1]=="Error"){
	#print(results)
	cat("An unexpected error occured during the simulation.  Ensure that the values of divide.data, pre.data and post.data provided are correct and re-run the functon.\n\n--- Simulation aborted ---\n\n")
	stop("An error occured during the simulation")
}
if(any(names(results)=="pilot.mcmc")){
	warning("The chains did not achieve convergence during the simulation, you should interpret the Bayesian results with extreme caution")
	converged <- FALSE
	results$mcmc <- results$pilot.mcmc
	results$summary <- results$pilot.summary
}else{
	converged <- TRUE
}

}else{
	results <- list(mcmc="MCMC analysis not performed")
}


blankquant <- quantile(0, probs=c(lci, 0.5, uci))

# No bootstrapping or WAAVP for paired model type data

if(any(c(dim(post)[1:2], dim(pre)[1:2])>1)){
	cat("Bootstrap and WAAVP calculations are not available for repeated pre and/or post treatment egg counts\n")
	boot.reductions = bootquant = method.boot = waavpquant = method.waavp =bootprob <- NA
}else{
	
	pre <- apply(pre.data, 3, mean) # Should only be 1 datapoint
	post <- apply(post.data, 3, mean) # Should only be 1 datapoint
	# Just to check:
	pre2 <- pre.data[1,1,]
	post2 <- post.data[1,1,]
	if(any(pre!=pre2) | any(post!=post2)) stop("An unexpected error occured while manipulating the data for the bootstrap/WAAVP methods")
	
	if(sum(control.animals)>0){
		warning("Control animals are ignored using the bootstrap and WAAVP calculations")
		pre <- pre[!control.animals]
		post <- post[!control.animals]
	}
	
	cat("Calculating the Bootstrap and WAAVP method analysis...\n")
	# Bootstrap:
	
	boot.pre <- matrix(data=sample(pre, length(pre)*bootstrap.iters, replace=TRUE), ncol=bootstrap.iters, nrow=length(pre))
	boot.post <- matrix(data=sample(post, length(post)*bootstrap.iters, replace=TRUE), ncol=bootstrap.iters, nrow=length(post))
	
	boot.reductions <- (1 - apply(boot.post, 2, mean) / apply(boot.pre, 2, mean)) *100

	boot.reductions[boot.reductions==-Inf] <- NA
	bootprob <- sum(boot.reductions < efficacy, na.rm=TRUE) / sum(!is.na(boot.reductions)) * 100
	
	bootquant <- quantile(boot.reductions, probs=c(lci, 0.5, uci), na.rm=TRUE)
	if(is.na(bootprob)){
		method.boot <- "Returned error"
	}else{
		method.boot <- if(bootprob >= (uci*100)) "Confirmed resistant" else if(bootprob >= 50) "Probable resistant" else if(bootprob >= (lci*100)) "Possible resistant" else "Confirmed susceptible"
	}

	# WAAVP:
	
	waavpquant <- blankquant
	
	pre.data <- pre
	post.data <- post

	N <- length(pre.data)
	arith.pre <- mean(pre.data)
	arith.post <- mean(post.data)
	var.pre <- var(pre.data)
	var.post <- var(post.data)

	red <- 100*(1-arith.post/arith.pre)
	var.red <- var.post/(N*arith.post^2) + var.pre/(N*arith.pre^2)
	try(if(var.post==0 & arith.post==0) var.red <- var.pre/(N*arith.pre^2))

	upper.ci <- 100 * (1-(arith.post/arith.pre * exp(-2.048*sqrt(var.red))))
	lower.ci <- 100 * (1-(arith.post/arith.pre * exp(2.048*sqrt(var.red))))

	waavpquant[1] <- lower.ci
	waavpquant[2] <- red
	waavpquant[3] <- upper.ci
	
	if(any(is.na(c(red,efficacy,lower.ci)))){
		method.waavp <- "Returned error"
	}else{
	method.waavp <- "susceptible"
	if(red < 95 | lower.ci < 90){
		method.waavp <- "Suspected resistant"
		if(red < 95 & lower.ci < 90){
			method.waavp <- "Confirmed resistant"
		}
	}
	
	if(confidence!=95){
		waavpquant[] <- NA
		method.waavp <- "Not available"
	}
}
}

if(!skip.mcmc){

cat("Calculating the Bayesian method analysis...\n")

reduction <- summary(results, vars='reduction(%)')
pre.mean <- summary(results, vars='groupmean')

mcred <- as.numeric(combine.mcmc(results, vars='reduction(%)'))
mcmcprob <- sum(mcred < efficacy) / length(mcred) * 100

mcmcquant <- blankquant
mcmcquant[2] <- reduction[,'Mean']
mcmcquant[1] <- reduction[,paste('Lower',confidence*100, sep='')]
mcmcquant[3] <- reduction[,paste('Upper',confidence*100, sep='')]

meanquant <- blankquant
meanquant[2] <- pre.mean[,'Mean']*divide.data
meanquant[1] <- pre.mean[,paste('Lower',confidence*100, sep='')]*divide.data
meanquant[3] <- pre.mean[,paste('Upper',confidence*100, sep='')]*divide.data

mcmcprob <- sum(reduction < efficacy) / length(reduction) * 100

method.mcmc <- if(mcmcprob >= (uci*100)) "Confirmed resistant" else if(mcmcprob >= 50) "Probable resistant" else if(mcmcprob >= (lci*100)) "Possible resistant" else "Confirmed susceptible"

#class <- (prob.res > 0.025) + (prob.res > 0.50) + (prob.res > 0.975) + 1
#method3 <- switch(class, "1"="susceptible", "2"="possible", "3"="probable", "4"="resistant")

}else{
	method.mcmc = mcmcquant = mcmcprob = meanquant = converged = efficacy <- NA
}

cat("Finished calculations\n")

output <- list(results.mcmc=method.mcmc, quant.mcmc=mcmcquant, results.boot=method.boot, quant.boot=bootquant, results.waavp=method.waavp, quant.waavp = waavpquant, prob.mcmc=mcmcprob, prob.boot=bootprob, meanquant=meanquant, confidence=confidence, converged=converged, efficacy=efficacy, animal.names=names, name=name)
if(record.chains) output <- c(output, list(mcmc=coda::as.mcmc.list(results)))

#browser()

if(write.file){
	cat("Writing results to file\n")
	filename <- new_unique(paste(name, ".results",sep=""), ".txt", ask=arguments$interactive)
	print.fecrt.results(output, filename=filename)
	filename <- new_unique(paste(name, ".graph",sep=""), ".pdf", ask=arguments$interactive)
	pdf(file=filename)
	if(!skip.mcmc) hist(mcred, col="red", breaks="fd", main="Posterior distribution for the true reduction", freq=FALSE, ylab="Probability", xlab="True FEC reduction (%)", xlim=c(0, 100))
	abline(v=efficacy)
	dev.off()
	if(record.chains) save(results, file=new_unique(paste("fecrt.", name, sep=""), ".Rsave", ask=FALSE))
	cat("Results file written successfully\n")
}else{
	if(plot.graph){
		if(!skip.mcmc){
			hist(mcred, col="red", breaks="fd", main="Posterior distribution for the true reduction", freq=FALSE, ylab="Probability", xlab="True FEC reduction (%)", xlim=c(0, 100))
			abline(v=efficacy)
		}
	}
}



cat("\nAnalysis complete\n\n", sep="")
cat("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*\n\n")
cat("--- End ---\n\n")

class(output) <- "fecrt.results"
return(output)


}



FECRT.analysis <- fecrt.analysis
fecrt <- fecrt.analysis
FECRT <- fecrt.analysis

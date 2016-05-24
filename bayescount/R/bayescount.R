bayescount <- function(name = NA, data = NA, setnames = NA, model = c("ZILP"), divide.data = 1, scale.mean = divide.data, remove.all.zeros = TRUE, test = TRUE, alt.prior = FALSE, write.file = TRUE, adjust.zi.mean = FALSE, likelihood=FALSE, record.chains=FALSE, ...)
{

	warning('The bayescount function is deprecated and will be removed from version 1.0 of bayescount (expected to be released in mid 2015)')

datanames <- setnames
omit.zeros <- remove.all.zeros
runname <- name

div <- divide.data  # Need both to prevent list data being divided twice

testwritable <- new_unique("test")
if(testwritable=="Directory not writable"){
	cat("\nThe working directory is not writable.  Please change the working directory\n\n")
	stop("Directory not writable")
}


passthrough <- list(...)

if(!is.null(passthrough$raw.output)) stop("Unable to return raw output of models when using the bayescount function.  Either use the bayescount.single function or use the record.chains option to write raw output to file for each dataset")

if(is.null(passthrough$silent.jags)) silent.jags <- eval(formals(autorun.jags)$silent.jags) else silent.jags <- passthrough$silent.jags
if(is.null(passthrough$jags)) jags <- eval(formals(autorun.jags)$jags) else jags <- passthrough$jags

test.jags <- testjags(jags, silent=TRUE)
if(test.jags[[2]][1]==FALSE){
	cat("Unable to call JAGS using '", jags, "' - try specifying the path to the JAGS binary as the jags argument\n", sep="")
	stop("Unable to call JAGS")
}

#now moved to autorun.jags:
#if((crash.retry != as.integer(crash.retry)) | (crash.retry < 0)){
#	cat("\nThe value of 'crash.retry' is not valid.  Please provide a non-negative integer\n\n")
#	stop("Parameter invalid")
#}

model <- toupper(model)

for(i in 1:length(model)){
	model[i] <- switch(model[i], P="SP", ZIP="ZISP", model[i])
}
models <- c("SP", "ZISP", "GP", "ZIGP", "WP", "ZIWP", "LP", "ZILP", "IP")
modelsfull <- c("single Poisson", "zero-inflated single Poisson", "gamma Poisson", "zero-inflated gamma Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson", "independant Poisson")

if(class(alt.prior)=="character"){
	stop("'alt.prior' must be either TRUE or FALSE for bayescount().  Use bayescount.single() for custom priors")
}


remove.missing <- TRUE # Partly removed facility to remove missing data since it doesn't make sense, and run.model removes all NA datasets now


if(model[1]=="ALL") model <- models

if((length(model) == 0) | (length(model) > length(models)) | sum(is.na(model)) > 0){
	cat("Invalid model selection.  Please choose from the following models: ", sep="")
	cat(models, sep=", ")
	cat("\n")
	stop("Invalid model selection")
}

for(i in 1:length(model)){
	if(sum(model[i]==models) != 1){
		cat("Invalid model selection \'", model[i], "\'.  Please choose from the following models: ", sep="")
		cat(models, sep=", ")
		cat("\n")
		stop("Invalid model selection")
	}
}

duplicated <- ""
for(i in 1:length(models)){
	if(sum(models[i]==model) > 1){
		cat("Duplicated model selection \'", models[i], "\'.  The model will be run only once!\n", sep="")
		duplicated <- c(duplicated, models[i])
	}
}

if(length(duplicated) > 1){
	for(i in 2:length(duplicated)){
		model[model==duplicated[i]] <- ""
		model <- c(model, duplicated[i])
	}
	model <- model[model!=""]	
}

models.nos <- 1:length(models)

model.no <- numeric(length=length(model))
modelfull <- character(length=length(model))

for(i in 1:length(model)){
	model.no[i] <- models.nos[models[models.nos]==model[i]]
	modelfull[i] <- modelsfull[model.no[i]]
}

if(sum(model=="ZISP") > 0 | sum(model=="ZIGP") > 0 | sum(model=="ZILP") > 0 | sum(model=="ZIWP") > 0 | sum(model=="ZIIP") > 0){
	any.zero.inflation <- TRUE
}else{
	any.zero.inflation <- FALSE
}

if(class(data)=="function") stop("The class of the object supplied for data is a function")



if(class(data)=="list"){
	
	if(suppressWarnings(class(data$totals) != "NULL" & class(data$repeats) != "NULL")){
		# data is a list of counts and repeats
		totals <- as.matrix(data$totals / divide.data)
		repeats <- as.matrix(data$repeats)
		
		div <- 1 # This stops the list type data being divided twice
		
		if(length(totals[1,])!=length(repeats[1,])) stop("The number of datasets is not equal between repeats and totals")
		if(length(totals[,1])!=length(repeats[,1])) stop("The number of repeats does not match the length of totals")
		
		n.datasets <- length(totals[1,])
		
		data <- array(NA, dim=c(length(totals[,1]), max(repeats), n.datasets))
		
		for(j in 1:n.datasets){
		for(i in 1:length(totals[,j])){
			if(is.na(totals[i,j]) & is.na(repeats[i,j])) next
			if(round(totals[i,j])!=totals[i,j]) stop("All totals must be an integer (after division by divide.data)")
			if(round(repeats[i,j])!=repeats[i,j] | repeats[i,j] < 1) stop("All repeats must be an integer greater than 0")
			
			if(repeats[i,j] > 1){

				done <- FALSE
	
				while(!done){
					data[i,1:repeats[i,j],j] <- rpois(repeats[i,j], totals[i,j]/repeats[i,j])
					if(sum(data[i,,j],na.rm=TRUE)==totals[i,j] & all(na.omit(data[i,,j]) >= 0)) done <- TRUE
		
				}
			}else{
				data[i,1,j] <- totals[i,j]
			}
		}
		}
		
	}else{
		
		stop("Data was provided as a list without elements 'totals' and 'repeats'.  To use a list of repeat samples from a single dataset, use bayescount.single()")
	}
}
	
if(class(data)=="array" & length(dim(data)) != 3) stop(paste("Data was provided as an array with ", length(dim(data)), " dimensions.  Arrays must have 3 dimensions (even if the third dimension is of length 1", sep=""))

if(class(data)=="numeric" | class(data)=="integer") data <- as.matrix(data)

cat("\n--- 'Bayescount': Analyse count data using a Bayesian distributional simulation model implemented in JAGS ---\n\n")
cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n")
cat("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*\n\n")

if(is.na(runname)){
	runname <- ask(prompt = "Please enter a name for this analysis:  ", type="character")
}

if(class(data)=="array" & identical(setnames[1],TRUE)) stop("Setnames cannot be TRUE for data provided as a list or array; please provide setnames manually or allow them to be generated automatically")

datain <- data
datana <- data
datana <- as.vector(na.omit(datana))
length(datana) <- 1

dataok=FALSE

setnamesheader <- identical(setnames, TRUE)

if(class(data)=="array" | class(data)=="matrix"){
	dataok <- TRUE
}else{
	if(is.na(datana)==FALSE){
	exists <- try(file.exists(datana), silent=TRUE)
	if((class(exists)=="try-error")==FALSE){
		if(exists==TRUE){
			suppressWarnings(data <- try(read.csv(datain, header=setnamesheader), silent=TRUE))
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
		data <- try(read.csv(datain, header=FALSE), silent=TRUE)
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



namesdone <- FALSE

if(identical(setnames[1],TRUE)){
	if(is.null(dimnames(data))){
		setnames <- data[1,]
		data <- matrix(data[2:nrow(data),], ncol=ncol(data))
	}else{
		setnames <- dimnames(data)[[length(dim(data))]]
	}
	names <- setnames
	namesdone <- TRUE
}else{
	if(identical(setnames[1],FALSE)){
		names <- paste("Dataset ", 1:dim(data)[length(dim(data))], sep="")
		namesdone <- TRUE
	}
	if(is.na(setnames[1])){
		if(class(dimnames(data)[[length(dim(data))]]) == "character"){
			names <- dimnames(data)[[length(dim(data))]]
		}else{
			names <- paste("Dataset ", 1:dim(data)[length(dim(data))], sep="")
		}
		namesdone <- TRUE
	}
}

if(!namesdone){
	if(length(setnames) != dim(data)[length(dim(data))]) stop("The length of the character vector setnames does not match the lenth of the data provided")
	names <- setnames
}

if(class(data)=="data.frame") data <- matrix(unlist(data), nrow=dim(data)[1])

if(class(data) == "matrix")	data <- array(data, dim=c(nrow(data), 1, ncol(data)), dimnames=list(NULL, NULL, dimnames(data)[[2]]))


if(class(data)!="array" | length(dim(data)) != 3) stop("An error occured while transforming the data to an appropriate format")

rownames <- FALSE # no longer using rownames



if(FALSE){  # ALL NOT BEING USED

names <- NA

if(namesdone==FALSE){
	if(length(setnames) != 1 | setnames[1] != FALSE){
		if(length(setnames) != dim(data)[3]) stop("The length of the character vector setnames does not match the lenth of the data provided")
		names <- setnames
		namesdone <- TRUE
	}else{
		if(class(dimnames(data)[[3]]) == "character"){
			names <- dimnames(data)[[3]]
		}else{
			names <- paste("Dataset ", 1:dim(data)[[3]], sep="")
		}
	}
}



if(is.na(datanames[1])==FALSE){
	if(datanames[1]!=TRUE){
		if(datanames[1]!=FALSE){
			if(rownames==FALSE){
				if(length(datanames)!=dim(data)[3]){
					stop("The length of the character vector setnames does not match the lenth of the data provided")
				}
			}else{
				if(length(datanames)!=(dim(data)[3]-1)){
					stop("The length of the character vector setnames does not match the lenth of the data provided")
				}
			}
			names <- datanames
			datanames <- TRUE
		}
	}
}

if(is.na(datanames[1])==TRUE){
	suppressWarnings(try({
		if(any(is.na(as.numeric(na.omit(data[1,1,]))))){
			names <- data[1,1,]
			data <- array(data[2:length(data[,1,1]),,], dim=c(dim(data)[1]-1, dim(data)[2], dim(data)[3]))
			datanames <- TRUE
		}else{
			if(class(dimnames(data)[[3]]) == "character"){
				names <- dimnames(data)[[3]]
				datanames <- TRUE
			}else{
				datanames <- FALSE
			}
		}
	}, silent=TRUE))
}

if(is.na(names)) names <- paste("Dataset ", 1:dim(data)[[3]], sep="")
if(FALSE){  # Not using this since rownames=FALSE and data is now an array
if(is.na(names[1])==FALSE){
	if(rownames==TRUE){
		ind.names <- data[ ,1]
		data <- data[ , 2:length(data[1,])]
	}else{
		ind.names <- paste("Row ", 1:length(data[,1]), sep="")
	}
}else{
	if(rownames==TRUE){
		if(datanames==TRUE){
			ind.names <- data[2:length(data[,1]),1]
			names <- data[1,2:length(data[1,])]
			data <- data[2:length(data[,1]), 2:length(data[1,])]
		}else{
			ind.names <- data[ ,1]
			names <- 1:(length(data[1,])-1)
			names <- paste("Dataset ", names, sep="")
			data <- data[ , 2:length(data[1,])]
		}
	}else{
		if(datanames==TRUE){
			ind.names <- paste("Row ", 1:(length(data[,1])-1), sep="")
			names <- data[1,]
			data <- data[2:length(data[,1]), ]
		}else{
			ind.names <- paste("Row ", 1:length(data[,1]), sep="")
			names <- 1:length(data[1,])
			names <- paste("Dataset ", names, sep="")
		}
	}

}
}
}  # END NOT BEING USED


dimensions <- dim(data)
suppressWarnings(data <- as.numeric(data) / div)
data <- array(as.numeric(data), dim=dimensions)

animalsindata1 <- length(data[,1,1])

names <- as.matrix(names)

if(div==1 & length(data) > 1){
	if(all(round(data/2)==(data/2)) | all(round(data/3)==(data/3))){
		cat("Warning:  It appears as though the data provided have been multiplied by a constant\n")
		if(test==TRUE){
			if(ask(prompt="Exit program to change 'divide.data'?", type="logical")) stop("User exited the program")
		}
	}
}
			
#print(names)
#print(data)
#stop()

if(rownames==TRUE && remove.missing==TRUE){  # Also not being used at the moment (rownames=FALSE)
	cat("Missing data cannot be removed if rownames are used.  Missing data will not be removed\n")
	remove.missing <- FALSE
}

cat("Settings are as follows:")
cat(paste("\nName of analysis:                             ", runname, sep=""))
#cat(paste("\nFirst row used for dataset names:             ", datanames, sep=""))
cat(paste("\nFirst column used for individual names:       ", rownames, sep=""))
cat(paste("\nNumber of datasets:                           ", as.numeric(length(names)), sep=""))
cat(paste("\nNumber of animals in dataset one              ", as.numeric(animalsindata1), sep=""))
cat(paste("\nDivide data by                                ", divide.data, sep=""))
cat(paste("\nAdjust mean by                                ", scale.mean, sep=""))

if(length(model)==1){
	cat(paste("\nModel to use:                                 ", modelfull[1], sep=""))
}else{
	cat(paste("\nModels to use:                                ", modelfull[1], " model", sep=""))
	for(i in 2:length(model)){
		cat(paste("\n                                              ", modelfull[i], " model", sep=""))
	}
}
cat(paste("\nOmit datasets with all zero counts:           ", omit.zeros, sep=""))
#cat(paste("\nRemove missing data:                          ", remove.missing, sep=""))
#cat(paste("\nMaximum time to spend per dataset:            ", max.time, sep=""))
#cat(paste("\nInteractive mode for autorun.jags:            ", interactive, sep=""))
#cat(paste("\nInitial burnin period:                        ", startsample, sep=""))

suppressWarnings(if(sum(model==c("GP", "ZIGP", "LP",	 "ZILP", "WP", "ZIWP"))>0){
	cat(paste("\nUse alternative prior for variance :          ", alt.prior, sep=""))
})
if(any.zero.inflation==TRUE){
	cat(paste("\nAdjust z-i means to mean of whole population: ", adjust.zi.mean, sep=""))
}
cat(paste("\nCalculate the likelihood for each model:      ", likelihood, sep=""))
cat(paste("\nSuppress JAGS output to screen:               ", silent.jags, sep=""))
#cat(paste("\nNumber of times to retry after crash:         ", as.numeric(crash.retry), sep=""))
cat(paste("\nWrite the results to file:                    ", write.file, sep=""))
cat(paste("\nSystem call to activate JAGS:                 ", jags, "\n\n", sep=""))
test.jags <- testjags(jags, silent=FALSE)


if(test==TRUE){
proceed <- FALSE
setdata <- (data[,,1])
cat("\n\nTesting the model function\n")
s1 <- try(testmod <- run.model(model = model[1], jags = jags, data = setdata, alt.prior = alt.prior, silent.jags = TRUE, summarise=FALSE))
s2 <- try(output <- extend.jags(testmod, burnin=1000, sample=1000, jags = jags, silent.jags = TRUE, summarise=FALSE))

if(class(s1)=="try-error" || class(s2)=='try-error'){
	proceed <- ask(prompt = "The test (first) dataset returned an error.  Continue with other datasets?", type="logical")
}else{
	cat("Test completed successfully\n")
	proceed <- TRUE
}
if(!proceed) stop("The program was halted by the user")
}

all.models <- model
all.modelfull <- modelfull
GP.largeod = ZIGP.largeod = WP.largeod = ZIWP.largeod = LP.largeod = ZILP.largeod <- 0

start.time <- Sys.time()

total.crashed <- 0
total.unconverged <- 0
total.error <- 0
total.largeod <- 0

cat("\n")

for(j in 1:length(all.models)){

	model <- all.models[j]
	modelfull <- all.modelfull[j]

	headers <- c("converged", "error.code", "samples", "samples.to.conv")
	resultscol <- 4

	headers <- c(headers, "mean.l.95", "mean.median", "mean.u.95")
	resultscol <- resultscol + 3
	
	headers <- c(headers, "coeff.variation.l.95", "coeff.variation.median", "coeff.variation.u.95")
	resultscol <- resultscol + 3
	
	if(model=="ZIGP" | model=="ZILP" | model=="ZIWP" | model=="ZISP" | any(strsplit(model, "")[[1]][2]==1:10)){
		headers <- c(headers, "zi.l.95", "zi.median", "zi.u.95")
		resultscol <- resultscol + 3
	}
	
	if(model=="GP" | model=="ZIGP" | model=="WP" | model=="ZIWP" | strsplit(model, "")[[1]][1]=="G"){
		headers <- c(headers, "scale.l.95", "scale.median", "scale.u.95")
		resultscol <- resultscol + 3
		headers <- c(headers, "shape.l.95", "shape.median", "shape.u.95")
		resultscol <- resultscol + 3
	}
	
	if(model=="LP" | model=="ZILP" | strsplit(model, "")[[1]][1]=="L"){
		headers <- c(headers, "log.mean.l.95", "log.mean.median", "log.mean.u.95")
		resultscol <- resultscol + 3
		headers <- c(headers, "log.variance.l.95", "log.variance.median", "log.variance.u.95")
		resultscol <- resultscol + 3
	}
	
	headers <- c(headers, "multivariate.psrf")
	resultscol <- resultscol + 1
	
	if(likelihood==TRUE){
		headers <- c(headers, "likelihood.l.95", "likelihood.median", "likelihood.u.95", "likelihood.MAX.OBS")
		resultscol <- resultscol + 4
	}
	
	headers <- c(headers, "time.taken")
	resultscol <- resultscol + 1
	
	zeros <- character(length=resultscol)
	zeros[] <- NA
	zeros[2] <- 3


	results <- matrix(NA, ncol=resultscol, nrow=length(names), dimnames =list(names, headers))
	headers <- paste(paste(headers, collapse=","), "\n", sep="")
	crashed <- 0
	unconverged <- 0
	error <- 0
	largeod <- 0
	
	name <- new_unique(name=paste(runname, ".", model, sep=""), suffix=".csv", ask = test*write.file, prompt="A results file with this name already exists.  Overwrite?")
	outfile <- file(name, 'w')
	cat("dataset,", headers, file = outfile, sep = "")
	
	for(i in 1:length(names)){
		done <- FALSE
		setdata <- data[,,i]
		
		if(remove.missing==TRUE){
			#setdata <- as.matrix(na.omit(data[,,i]))  Handled by run.model now
		}else{
			#setdata <- as.integer(data[,,i])
		}
		if(omit.zeros==TRUE){
			if(sum(setdata,na.rm=TRUE) < 1){
				results[i,] <- zeros
				output <- zeros
				cat("Dataset '", names[i], "' contained all zeros and was therefore omitted\n", sep="")
				done <- TRUE
			}
		}
		tries.left <- 1#crash.retry + 1 #crash.retry moved to autorun.jags
		pre.time <- Sys.time()
		while(done==FALSE){
			cat("\nRunning the ", model, " model for dataset '", names[i], "'...\n", sep="")
			#if(length(setdata)!=length(data[,,i]) && remove.missing==TRUE){
			#	cat("\n*WARNING*  ", length(data[,1,i]) - length(setdata[,1]), " missing and/or non-numeric datapoints were removed from dataset '", names[i], "'\n", sep="")
			#}    Missing data removed by run.model now
			#if(sum(is.na(setdata[,1])) > 0 && remove.missing==FALSE){
			#	cat("\n*WARNING*  There are ", sum(is.na(setdata)), " missing and/or non-numeric datapoints in dataset '", names[i], "'\n", sep="")
			#}
			output <- bayescount.single(model=model, data = setdata, alt.prior = alt.prior, adjust.zi.mean = adjust.zi.mean, likelihood=likelihood, raw.output=if(record.chains){list(name=paste(strsplit(name, ".csv")[[1]], collapse=".csv"), setname=names[i])}else{FALSE}, silent.jags=silent.jags, ...)
			
			if(((output[2]=="1") | (output[2]=="5")) && (sum(na.omit(output[1] == "1")) == 0)){
				cat("The ", model, " model for dataset '", names[i], "' returned an error", sep="")
				tries.left <- tries.left-1
				if(tries.left==0){
					cat("\n\n")
					done <- TRUE
				}else{
					cat(".  Re-trying:\n") # no longer used
				}
			}else{
				done <- TRUE
			}
		}
		post.time <- Sys.time()
		
#		if(output[2]=="1"){
#			#crashed <- crashed + 1 #Getting rid of crashed for now
#			#total.crashed <- total.crashed + 1
#			error <- error + 1
#			total.error <- total.error + 1
#		}
		if(output[2]=="1" | output[2]=="5"){
			error <- error + 1
			total.error <- total.error + 1
		}
		if(sum(na.omit(output[1]=="0")) == 1 && (sum(na.omit(setdata)) > 0 || omit.zeros==FALSE)){
			if(output[2]=="2"){
				cat("*WARNING*  The ", model, " model for dataset '", names[i], "' failed to achieve convergence\n", sep="")
				unconverged <- unconverged + 1
				total.unconverged <- total.unconverged + 1
			}else{
				cat("*WARNING*  The PSRF for the ", model, " model for dataset '", names[i], "' increased above the convergence target for the final chains\n", sep="")
			}
		}
		
		if(sum(na.omit(output[1] == "1")) == 1){			
			largeod <- largeod + as.numeric(assess.variance(model=list(model=model, alt.prior=alt.prior, l.95=output["coeff.variation.l.95"], u.95=output["coeff.variation.u.95"])))
		}

		success <- try(results[i,] <- as.numeric(output))
		if(class(success)=="try-error"){
			if(!any(output[2]==1:5)){ # pick up any errors that are missed because of strange returns
				error <- error + 1
				total.error <- total.error + 1
			}
			
			cat("ERROR: Expecting ", length(results[i,]), " values representing ", headers, " and got '", paste(names(output), collapse=", "), "'.\n", sep="")
			warning(paste("Expecting ", length(results[i,]), " values representing ", headers, " and got '", paste(names(output), collapse=", "), "' for dataset ", names[i], " with the ", model, " model.", sep=""))
		}else{
			if(omit.zeros==FALSE | sum(setdata,na.rm=TRUE) > 0){
				results[i,c("mean.l.95", "mean.median", "mean.u.95")] <- results[i,c("mean.l.95", "mean.median", "mean.u.95")] * scale.mean
			}
		}
		time.taken <- timestring(pre.time, post.time, show.units=TRUE)
		total.time <- timestring(start.time, post.time, show.units=TRUE)
		number.models <- length(all.models) * length(names)
		current.number <- ((j-1) * length(names)) + i
		percent.complete <- current.number / number.models
		time.remaining <- timestring((as.integer(difftime(post.time, start.time, units="secs")) / percent.complete) - as.integer(difftime(post.time, start.time, units="secs")), show.units=TRUE)
		
		if(sum(na.omit(setdata)) > 0 || omit.zeros==FALSE){
			cat("Dataset '", names[i], "' completed in ", time.taken, " with the ", model, " model\n", sep="")
		}
		cat(i+(length(names)*(j-1)), " of ", number.models, " datasets (", round(percent.complete, digits=2)*100, "%) completed\n", sep="")
		cat(total.time, " elapsed, estimated time remaining: ", time.remaining, "\n", sep="") 
		if(j > 1){
			cat(unconverged, " failed convergence, and ", error, " quit with an error with the ", model, " model\n", sep="")
			cat(total.unconverged, " failed convergence, and ", total.error, " quit with an error for all models\n", sep="")
		}else{
			cat(unconverged, " failed convergence, and ", error, " quit with an error\n", sep="")
		}
		cat("\n")
		cat(names[i], results[i,], file = outfile, sep = ",")
		cat("\n", file=outfile, sep="")
		gc()
	}
	close(outfile)
	if(write.file==FALSE){
		unlink(name)
	}
#	assign(paste(runname, ".", model, ".results", sep=""), results, pos=".GlobalEnv")
	if((model=="GP") | (model=="ZIGP") | (model=="WP") | (model=="ZIWP")){
		assign(paste(model, ".largeod", sep=""), largeod)
	}
}


cat("All models completed.  Total time taken:  ", timestring(start.time, post.time), "\n\n", sep="")
cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n")
cat("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*\n\n")
cat("--- End ---\n\n")


largeods <- c(GP.largeod, ZIGP.largeod, WP.largeod, ZIWP.largeod, LP.largeod, ZILP.largeod)
odmodels <- c("gamma Poisson", "zero-inflated gamma Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson")

for(i in 1:length(largeods)){
	if((alt.prior==TRUE | is.character(alt.prior)) && largeods[i]>0){
		cat("*WARNING*  The 95% confidence interval for the variance was very large in a total of ", largeods[i], " datasets using the ", odmodels[i], " model, indicating a lack of information about this parameter in these datasets.  These models should be re-run with the standard prior distribution and the two sets of results compared to ensure that the results for all parameters are reliable\n")
	}
	if(alt.prior==FALSE && largeods[i]>0){
		cat("*WARNING*  The 95% confidence interval for the variance was very large in a total of ", largeods[i], " datasets using the ", odmodels[i], " model, indicating a lack of information about this parameter in these datasets.  These models should be re-run with an alternative prior distribution and the two sets of results compared to ensure that the results for all parameters are reliable\n")
	}
}

}
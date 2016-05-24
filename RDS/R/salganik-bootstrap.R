salganik.bootstrap.estimates <- function(rds.data, group.variable,number.of.bootstrap.samples,
		estimator.name, N=NULL,...){
	
	if(is(rds.data,"rds.data.frame")){
		stopifnot(group.variable %in% names(rds.data))
		network.size <- attr(rds.data, "network.size.variable")
	}else{
		stop("rds.data must be of type rds.data.frame")
	}
	
	############################################################################################
	# First extract the needed information from the rds data frame object. 
	group.memberships <- as.vector(rds.data[[group.variable]])
	group.names <- unique(group.memberships)
	
	# NB: if one wants to treat missing values as a separate group,
	# then the following will need to be changed. 
	group.names <- group.names[!is.na(group.names)]
	
	number.of.groups <- length(group.names)
	
	# Now translate these to indices.  This is done for the sake of efficiency.
	group.indices <- match(group.memberships,group.names)
	
	# We also want to extract the wave information.  We will do this in order to 
	# sort the observations by wave.
	wave <- get.wave(rds.data)
	
	# If everything is ordered by ascending wave we can
	# simulate by running through the data set by rows.
	wave.order <- order(wave)    
	group.indices <- group.indices[wave.order]
	id <- get.id(rds.data)[wave.order]
	recruiter.id <- get.rid(rds.data)[wave.order]
	degrees <- data.frame(rds.data)[wave.order,network.size]    
	wave <- wave[wave.order]
	
	# Now rationalize the recruiter.id information so that it corresponds to recruiter
	# row.
	recruiter.row <- match(recruiter.id,id)
	
	# The zeros in the recruiter id will be mapped to NA's
	recruiter.row[is.na(recruiter.row)] <- 0
	
	wave.one.start <- match(1,wave)
	seed.rows <- 1:(wave.one.start - 1)
	number.of.seeds <- length(seed.rows)
	sample.size <- length(id)
	
	
	# Now we need a count of transitions.
	tij <- count.transitions(rds.data,group.variable)
	
	#if no observed transition, use the marginal
	cj <- colSums(tij)
	no.trans <- rowSums(tij)<.5
	tij[no.trans,] <- cj
	
	sample.size <- sum(!is.na(group.indices[wave > 0]))
	
	degree.i <- lapply(1:number.of.groups,
			function(g) {
				d <- degrees[group.indices == g]
				d[!is.na(d)]
			})
	################################################################################
	# This internal function creates a single bootstrap sample.  The return value 
	# is a data frame with id, recruiter.id, group variable and network size.  
	bootstrapper <- function(){
		seed.group.index <- sample(group.indices[!is.na(group.indices)],size=1)
		while(all(tij[seed.group.index,]==0)){
			seed.group.index <- sample(group.indices[!is.na(group.indices)],size=1)
		}
		
		results <- matrix(nrow=sample.size + 1, ncol=2)
		current.idx <- 1
		resample <- function(x, ...) {
			if(is.null(x))
				NA 
			else 
				x[sample.int(length(x), ...)]
		}
		results[1,1] <- seed.group.index
		results[1,2] <- resample(degree.i[[results[1,1]]],size=1)
		for(i in 2:(sample.size+1)){
			previous.group <- results[i-1,1]
			results[i,1] <- sample.int(number.of.groups,size=1,prob=tij[previous.group,])
			results[i,2] <- resample(degree.i[[results[i,1]]],size=1)
		}
		
		
		colnames(results) <- c('group.index','degree')
		
		bootstrapped.data <- data.frame(id=1:(sample.size+1),
				recruiter.id = 0:sample.size,
				network.size.variable=results[,2])
		bootstrapped.data[,group.variable] <- group.names[results[,1]]
		bootstrapped.data <- as.rds.data.frame(bootstrapped.data,
				population.size=get.population.size(rds.data),
				check.valid=FALSE)
		
		return(bootstrapped.data)  
	}
	
	f <- function(){
		RDS.estimates.local(
			rds.data=bootstrapper(),
			outcome.variable=group.variable,
			weight.type=estimator.name,
			empir.lik=FALSE,
			N=N,
			...)@estimate
	}
	
	bs.results <- replicate(number.of.bootstrap.samples, f())
	
	value <- matrix(0,ncol=length(group.names),nrow=number.of.bootstrap.samples)
	colnames(value) <- group.names
	colnames(value)[colnames(value)=="NA.NA"] <- "NA"
	if(is.matrix(bs.results)){
		for(i in 1:nrow(bs.results)){
			value[,i] <-  unlist(bs.results[i,])
		}
	}else{
		if(is.list(bs.results)){
			for(i in 1:number.of.bootstrap.samples){
				value[i,names(bs.results[[i]])] <-  unlist(bs.results[[i]])
			}
		}else{
			value[,1] <- bs.results
			value[,2] <- bs.results
		}
	}

	
	return(value) 
	
}




salganik.bootstrap.se <- function(rds.data,group.variable,
		number.of.bootstrap.samples,estimator.name,N=NULL,...)
{
	result <- salganik.bootstrap.estimates(
			rds.data=rds.data,
			group.variable=group.variable,
			number.of.bootstrap.samples=number.of.bootstrap.samples,
			estimator.name=estimator.name,
			N=N,
			...)
	
	result <- result[apply(!is.nan(result),1,all),]
	if(nrow(result)>1){
		a=(sqrt(diag(stats::var(result))))	
	}else{
		a=(sqrt(stats::var(as.numeric(result))))	
	}
	attr(a,"bsresult") <- result
	return(a)
}

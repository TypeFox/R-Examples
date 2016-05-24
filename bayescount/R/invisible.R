print.fecrt.results <- function(x=list(), filename=FALSE, ...){
		
	if(filename!=FALSE) file=filename else file=""
	
	N <- length(x$animal.names)
	
	digits <- 1
	
	if(filename==FALSE) cat("\n", file=file,append=TRUE)
	
	cat(paste("Results of the feacal egg count reduction test analysis for '", x$name, "' with ", N, " animals:\n\n", sep=""), file=file,append=TRUE)
	
	if(all(is.na(x$results.mcmc))) cat("The Bayesian MCMC analysis method was not used\n\n", file=file,append=TRUE) else cat(paste("The Bayesian MCMC method concluded that these data have a ", round(x$prob.mcmc, digits), "% probability of coming from a resistant herd, which is defined as '", x$results.mcmc, "'.  The ", x$confidence*100, "% credible interval and median result for the mean efficacy according to the Bayesian method is:\nl", x$confidence*100, ": ", round(x$quant.mcmc[1], digits), "%    median: ", round(x$quant.mcmc[2], digits), "%    u", x$confidence*100, ": ", round(x$quant.mcmc[3], digits), "%\n\n", sep=""), file=file,append=TRUE)
	
	if(identical(NA, x$results.boot)){
		cat("Results from bootstrapping and the WAAVP method were not available\n\n", file=file,append=TRUE)
		
	}else{
		cat(paste("The bootstrap method concluded that these data have a ", round(x$prob.boot,digits), "% probability of coming from a resistant herd, which is defined as '", x$results.boot, "'.  The ", x$confidence*100, "% confidence interval and median result for the mean efficacy according to the bootstrap method is:\nl", x$confidence*100, ": ", round(x$quant.boot[1], digits), "%    median: ", round(x$quant.boot[2],digits), "%    u", x$confidence*100, ": ", round(x$quant.boot[3],digits), "%\n\n", sep=""), file=file,append=TRUE)
		if(x$confidence!=0.95){
			cat("Results from the WAAVP method are not available unless confidence is 95%\n\n")
		}else{
			cat(paste("The WAAVP method defined the data as '", x$results.boot, "'.\nThe 95% confidence interval and median result for the mean efficacy according to the WAAVP method is:\nl95: ", round(x$quant.waavp[1],digits), "%    median: ", round(x$quant.waavp[2],digits), "%    u95: ", round(x$quant.waavp[3],digits), "%\n\n", sep=""), file=file,append=TRUE)
		}
	}
	
	if(all(!is.na(x$results.mcmc))){
	cat("\nThe following additional results were obtained using the Bayesian MCMC method:\n\n", file=file,append=TRUE)
	cat(paste("The median and ", x$confidence*100, "% credible intervals for the true mean egg count before treatment:\nl", x$confidence*100, ": ", round(x$meanquant[1],digits+1), "    median: ", round(x$meanquant[2],digits+1), "    u", x$confidence*100, ": ", round(x$meanquant[3],digits+1), "\n\n", sep=""), file=file,append=TRUE)
	#if(!all(is.na(x$indredquant))){
	#	cat(paste("Individual animal probabilities of 'efficacy < ", x$efficacy, "%', along with median and ", x$confidence, "% credible intervals for mean egg count reduction:\n\n", sep=""), file=file,append=TRUE)
	#for(i in 1:N){
	#	cat(paste(x$animal.names[i], ":  prob reduction < ", x$efficacy, "%: ", round(x$ind.prob.inf[i],digits), "%    l", x$confidence, ": ", round(x$indredquant[i,1],digits), "%    median:", round(x$indredquant[i,2],digits), "%    u", x$confidence, ": ", round(x$indredquant[i,3],digits), "%\n", sep=""), file=file,append=TRUE)
	#}
	#cat("\n", file=file, append=TRUE)
	#}else{
	#	cat("Individual animal efficacies were not assessed\n\n", file=file, append=TRUE)
	#}
	
	if(!x$converged) cat("*WARNING* The chains did not achieve convergence during the simulation, you should interpret the Bayesian MCMC results with extreme caution\n\n")
	
}

	cat("*THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n", file=file,append=TRUE)
	
    cat("\n", file=file,append=TRUE)
	
}

assess.variance <- function(model){
	
	alt.prior <- model$alt.prior
	l.95 <- model$l.95
	u.95 <- model$u.95
	model <- toupper(model$model)

	largeod <- FALSE

	if(model=="GP" | model=="ZIGP"){				
		if(alt.prior==TRUE | is.character(alt.prior)){
			if((l.95 < 0.002) && (u.95 > 0.1)){
				largeod <- TRUE
			}
		}else{
			if((l.95 < 0.01) && (u.95 > 100)){
				largeod <- TRUE
			}
		}
	}
	
	if(model=="WP" | model=="ZIWP"){				
		if(alt.prior==TRUE | is.character(alt.prior)){
			if((l.95 < 0.0001) && (u.95 > 0.01)){
				largeod <- TRUE
			}
		}else{
			if((l.95 < 0.001) && (u.95 > 10)){
				largeod <- TRUE
			}
		}
	}

	if(model=="LP" | model=="ZILP"){				
		if(alt.prior==TRUE | is.character(alt.prior)){
			if((l.95 < 0.01) && (u.95 > 1000)){
				largeod <- TRUE
			}
		}else{
			if((l.95 < 0.001) && (u.95 > 10)){
				largeod <- TRUE
			}
		}
	}

	return(largeod)
}

binary.search <- function(f, target, interval=c(-Inf, Inf), start=NA, precision=5){
	
	min <- interval[1]
	max <- interval[2]
	
	if(is.na(start)){
		if(min == -Inf) smin <- -10 else smin <- min
		if(max == Inf) smax <- 10 else smax <- max
		start <- runif(1, smin, smax)
	}
	
	value <- start
	i <- 0
	while(max==Inf){
		new <- f(value)
		if(new > target) max <- value
		if(new < target) min <- value
		value <- value + (value^2)*10 + as.integer(value==0)
		i <- i+1
		if(i > 100 | new==Inf) return(list(value=NA, objective=NA, status=paste("Max not found (up to ", min, ")", sep="")))
	}
	
	value <- start
	i <- 0
	while(min==-Inf){
		new <- f(value)
		if(new > target) max <- value
		if(new < target) min <- value
		value <- value - ((value^2)*10 + as.integer(value==0))
		i <- i+1
		if(i > 100 | new==-Inf) return(list(value=NA, objective=NA, status=paste("Min not found (down to ", max, ")", sep="")))
	}
	
	if(f(min)>target) return(list(value=NA, objective=NA, status="Minimum value is below that specified"))
	if(f(max)<target) return(list(value=NA, objective=NA, status="Maximum value is above that specified"))
	
	
	values <- numeric(10)
	news <- numeric(10)

	done <- FALSE
	while(!done){
		values <- seq(from=min, to=max, length=10)
		for(i in 1:10){
			news[i] <- f(values[i])
		}
		if(any(news > target)) max <- min(values[which(news>target)])
		if(any(news < target)) min <- max(values[which(news<target)])
		
		if(any(round(news, precision)==round(target,precision))){
			errors <- (news-target)^2
			new <- mean(news[which(errors==min(errors))])
			value <- mean(values[which(errors==min(errors))])
			status <- "OK"
			done <- TRUE
		}
		if(max-min < (.Machine$double.eps)^0.25){
			#newsminrmse <- news[which((news-target)^2==min((news-target)^2))]
			valuesminrmse <- values[which((news-target)^2==min((news-target)^2))]
			
			value <- mean(valuesminrmse)
			status <- "Absolute value not possible but this is closest"
			done <- TRUE
		}
		
	}
	return(list(value=value, objective=f(value), status=status))
}

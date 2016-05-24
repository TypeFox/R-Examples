maximise.likelihood <- function(data=stop("Data must be specified"), model=stop("Please specify a distribution"), mean=NA, variance=NA, zi=NA, shape=NA, scale=NA, silent=FALSE){
	
	warning('The maximise.likelihood function is deprecated and will be removed from version 1.0 of bayescount (expected to be released in mid 2015)')
	
	model <- toupper(model)
	testdata <- data
	
	models <- c("P", "ZIP", "G", "ZIG", "L", "ZIL", "W", "ZIW", "GP", "ZIGP", "LP", "ZILP", "WP", "ZIWP")
		
	if((length(model) != 1) |  sum(is.na(model)) > 0 | sum(model==models)!=1){
		if(silent==FALSE){
			cat("Invalid model selection.  Please choose from ONE of the following distributions: ", sep="")
			cat(models, sep=", ")
			cat("\n")
		}
		stop("Invalid model selection")
	}
	
	if(silent==FALSE) cat("\n")
	
	if(sum(is.na(data)) > 0){
		if(silent==FALSE) cat("WARNING:  Missing data was removed before calculating the likelihood\n\n")
		data <- na.omit(data)
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
	
	#if(guitest$os=='windows' | (guitest$R.GUI == "AQUA" & guitest$R.package.type == "mac.binary")) eol <- "" else eol <- "                  \n"
	if(guitest$R.GUI == "AQUA" & guitest$R.package.type == "mac.binary"){
		clearline <- ""
		eol <- "\n"
	}else{
		clearline <- "\r                                                         \r"
		eol <- "\t\t"
	}
	
	if(silent==FALSE){
		cat("Maximising the likelihood for the '", model, "' model.  This may take some time\n")
		
		flush.console()
	}
	
	if(model=="P"){
		if(is.na(mean)) mean <- mean(data)
		
		if(silent==TRUE){
			f1 <- function(mean) return(max(-Inf, likelihood(model=model, data=data, mean=mean, silent=as.character(TRUE)), na.rm=TRUE))
		}else{
			cat("\n\tmean\n", eol, "\t", signif(mean+(mean/100000), 4), eol, sep="")
			f1 <- function(mean){
				cat(clearline, "\t", signif(mean+(mean/100000), 4), eol, sep="")
				flush.console()
				return(max(-Inf, likelihood(model=model, data=data, mean=mean, silent=as.character(TRUE)), na.rm=TRUE))
			}
		}
		maxim <- optimise(f1, mean, lower=0, upper=max(data)*10, maximum=TRUE)
		results <- c(mean=maxim$maximum)
	}

	if(model=="ZIP"){
		#if(is.na(mean)) mean <- mean(data[data!=0])
		#if(is.na(zi)) zi <- sum(data==0)/length(data)*100
		
		if(is.na(mean)) mean <- mean(data)
		if(is.na(zi)) zi <- 0
		
		if(silent==TRUE){
			f2 <- function(pars=c(NA, NA)){
				if(pars[2] < 0 | pars[2] > 100 | pars[1] < 0) return(-Inf)
				return(max(-Inf, likelihood(model=model, data=data, mean=pars[1], zi=pars[2], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}else{
			cat("\n\tmean\t\tzi\n", eol, "\t", signif(mean+(mean/100000), 4), "\t\t", signif(zi+(zi/100000), 4), eol, sep="")
			f2 <- function(pars=c(NA, NA)){
				if(pars[2] < 0 | pars[2] > 100 | pars[1] < 0) return(-Inf)
				cat(clearline, "\t", signif(pars[1]+(pars[1]/100000), 4), "\t\t", signif(pars[2]+(pars[2]/100000), 4), eol, sep="")
				flush.console()
				return(max(-Inf, likelihood(model=model, data=data, mean=pars[1], zi=pars[2], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}
		maxim <- optim(c(mean, zi), f2, control=list(fnscale=-1))
		results <- c(mean=maxim$par[1], zi=maxim$par[2])
	}
	
	if(model=="G" | model=="GP"){
		if(is.na(scale)) scale <- var(data)/mean(data)
		if(is.na(shape)) shape <- mean(data) / scale
		
		if(silent==TRUE){
			f3 <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0) return(-Inf)
				return(max(-Inf, likelihood(model=model, data=data, shape=pars[1], scale=pars[2], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}else{
			cat("\n\tshape\t\tscale\n", eol, "\t", signif(shape, 4), "\t\t", signif(scale, 4), eol, sep="")
			f3 <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), eol, sep="")
				flush.console()
				return(max(-Inf, likelihood(model=model, data=data, shape=pars[1], scale=pars[2], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}
		maxim <- optim(c(shape, scale), f3, control=list(fnscale=-1))
		results <- c(shape=maxim$par[1], scale=maxim$par[2])
	}
	
	if(model=="ZIG" | model=="ZIGP"){
		
		#if(is.na(scale)) scale <- var(data[data!=0])/mean(data[data!=0])
		#if(is.na(shape)) shape <- mean(data[data!=0]) / scale
		#if(is.na(zi)) zi <- sum(data==0)/length(data)*100
		
		if(is.na(scale)) scale <- var(data)/mean(data)
		if(is.na(shape)) shape <- mean(data) / scale
		if(is.na(zi)) zi <- 0
		
		if(silent==TRUE){
			f4 <- function(pars=c(NA, NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				return(max(-Inf, likelihood(model=model, data=data, shape=pars[1], scale=pars[2], zi=pars[3], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}else{
			cat("\n\tshape\t\tscale\t\tzi\n", eol, "\t", signif(shape, 4), "\t\t", signif(scale, 4), "\t\t", signif(zi, 4), eol, sep="")
			f4 <- function(pars=c(NA, NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), "\t\t", signif(pars[3], 4), eol, sep="")
				flush.console()
				return(max(-Inf, likelihood(model=model, data=data, shape=pars[1], scale=pars[2], zi=pars[3], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}
		maxim <- optim(c(shape, scale, zi), f4, control=list(fnscale=-1))
		results <- c(shape=maxim$par[1], scale=maxim$par[2], zi=maxim$par[3])
	}
	
	if(model=="W" | model=="WP"){
		if(is.na(scale)) scale <- 4
		if(is.na(shape)) shape <- mean(data)
		
		if(silent==TRUE){
			f5 <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0) return(-Inf)
				return(max(-Inf, likelihood(model=model, data=data, shape=pars[1], scale=pars[2], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}else{
			cat("\n\tshape\t\tscale\n", eol, "\t", signif(shape, 4), "\t\t", signif(scale, 4), eol, sep="")
			f5 <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), eol, sep="")
				flush.console()
				return(max(-Inf, likelihood(model=model, data=data, shape=pars[1], scale=pars[2], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}
		maxim <- optim(c(shape, scale), f5, control=list(fnscale=-1))
		results <- c(shape=maxim$par[1], scale=maxim$par[2])
	}
	
	if(model=="ZIW" | model=="ZIWP"){
		if(is.na(scale)) scale <- 4
		if(is.na(shape)) shape <- mean(data[data!=0])
		if(is.na(zi)) zi <- sum(data==0)/length(data)*100
		
		if(silent==TRUE){
			f6 <- function(pars=c(NA, NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				return(max(-Inf, likelihood(model=model, data=data, shape=pars[1], scale=pars[2], zi=pars[3], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}else{
			cat("\n\tshape\t\tscale\t\tzi\n", eol, "\t", signif(shape, 4), "\t\t", signif(scale, 4), "\t\t", signif(zi, 4), eol, sep="")
			f6 <- function(pars=c(NA, NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), "\t\t", signif(pars[3], 4), eol, sep="")
				flush.console()
				return(max(-Inf, likelihood(model=model, data=data, shape=pars[1], scale=pars[2], zi=pars[3], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}
		maxim <- optim(c(shape, scale, zi), f6, control=list(fnscale=-1))
		results <- c(shape=maxim$par[1], scale=maxim$par[2], zi=maxim$par[3])
	}
	
	if(model=="L" | model=="LP"){
		if(is.na(mean)) mean <- mean(data)
		if(is.na(variance)) variance <- var(data)
		
		if(silent==TRUE){
			f7 <- function(pars=c(NA, NA)){
				if(pars[2] < 0 | pars[1] < 0) return(-Inf)
				return(max(-Inf, likelihood(model=model, data=data, mean=pars[1], variance=pars[2], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}else{
			cat("\n\tmean\t\tvariance\n", eol, "\t", signif(mean, 4), "\t\t", signif(variance, 4), eol, sep="")
			f7 <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] < 0) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), eol, sep="")
				flush.console()
				return(max(-Inf, likelihood(model=model, data=data, mean=pars[1], variance=pars[2], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}
		maxim <- optim(c(mean, variance), f7, control=list(fnscale=-1))
		results <- c(mean=maxim$par[1], variance=maxim$par[2])
	}
	
	if(model=="ZIL" | model=="ZILP"){
		#if(is.na(mean)) mean <- mean(data[data!=0])
		#if(is.na(variance)) variance <- var(data[data!=0])
		#if(is.na(zi)) zi <- sum(data==0)/length(data)*100
		
		if(is.na(mean)) mean <- mean(data)
		if(is.na(variance)) variance <- var(data)
		if(is.na(zi)) zi <- 0
		
		if(silent==TRUE){
			f8 <- function(pars=c(NA, NA, NA)){
				if(pars[2] < 0 | pars[1] < 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				return(max(-Inf, likelihood(model=model, data=data, mean=pars[1], variance=pars[2], zi=pars[3], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}else{
			cat("\n\tmean\t\tvariance\t\tzi\n", eol, "\t", signif(mean, 4), "\t\t", signif(variance, 4), "\t\t", signif(zi, 4), eol, sep="")
			f8 <- function(pars=c(NA, NA, NA)){
				if(pars[2] < 0 | pars[1] < 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), "\t\t", signif(pars[3], 4), eol, sep="")
				flush.console()
				return(max(-Inf, likelihood(model=model, data=data, mean=pars[1], variance=pars[2], zi=pars[3], silent=as.character(TRUE)), na.rm=TRUE))
			}
		}
		maxim <- optim(c(mean, variance, zi), f8, control=list(fnscale=-1))
		results <- c(mean=maxim$par[1], variance=maxim$par[2], zi=maxim$par[3])
	}
	
	if(silent==FALSE){
		#if(guitest$R.GUI == "AQUA" & guitest$R.package.type == "mac.binary") cat("\n")
		cat("\n\nFinished maximising the likelihood\n\n")
	}
	results <- c(results, log.likelihood=maxim$value)
	return(results)
}
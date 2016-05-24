



#' Calculates estimates at each successive wave of the sampling process
#' @param rds.data An rds.data.frame
#' @param outcome.variable The outcome
#' @param est.func A function taking rds.data and outcome.variable as parameters and
#' 	returning an rds.weighted.estimate object
#' @param ... additional parameters for est.func
cumulative.estimate <- function(rds.data, outcome.variable, est.func=RDS.II.estimates, ...){

	wave <- get.wave(rds.data)
	max.wave <- max(wave)
	if(is.character(rds.data[[outcome.variable]])){
		rds.data[[outcome.variable]] <- as.factor(rds.data[[outcome.variable]])
	}
	if(all(is.na(rds.data[[outcome.variable]]))){
		o <- rds.data[[outcome.variable]]
		if(is.factor(o)){
			estimates <- matrix(NA,nrow=max.wave+1,ncol=length(levels(o)))
			colnames(estimates) <- levels(o)	
		}else{
			estimates <- matrix(NA,nrow=max.wave+1,ncol=1)
		}
		rownames(estimates) <- 0:max.wave	
		attr(estimates,"n") <- rep(0,max.wave+1)
		return(estimates)
	}
	ests <- list()
	nw <- rep(NA,max.wave+1)
	for(w in 0:max.wave){
		wtd.est <- est.func(rds.data=rds.data[wave<=w,],
				outcome.variable=outcome.variable,
				empir.lik=FALSE,
				...)
		ests[[w+1]] <- wtd.est@estimate
		isna <- is.na(rds.data[wave<=w,outcome.variable])
		poswt <- wtd.est@weights > .Machine$double.eps
		nw[w+1] <- sum(!isna[wtd.est@subset] & poswt)

	}
	
	levs <- names(ests[[max.wave+1]])
	if(!is.null(levs))
		nlevs <- length(levs)
	else
		nlevs <- 1
	estimates <- matrix(0,nrow=max.wave+1,ncol=nlevs)
	colnames(estimates) <- levs
	rownames(estimates) <- 0:max.wave
	for(i in 1:(max.wave+1)){
		if(!is.null(levs)){
			inds <- match(names(ests[[i]]),levs)
			estimates[i,inds] <- ests[[i]]
		}else{
			estimates[i,1] <- ests[[i]]
		}
	}
	attr(estimates,"n") <- nw
	estimates
}



#' Convergence Plots
#' @description This function creates diagnostic convergence plots for RDS estimators.
#' @param rds.data An rds.data.frame.
#' @param outcome.variable A character vector of outcome variables.
#' @param est.func A function taking rds.data and outcome.variable as parameters and
#' 	returning an rds.weighted.estimate object.
#' @param ... additional parameters for est.func.
#' @param as.factor Convert all outcome variables to factors
#' @references 
#' Krista J. Gile, Lisa G. Johnston, Matthew J. Salganik \emph{Diagnostics for Respondent-driven Sampling} eprint arXiv:1209.6254, 2012 
#' @examples 
#' data(faux)
#' convergence.plot(faux,c("X","Y"))
#' @export
convergence.plot <- function(rds.data, outcome.variable, est.func=RDS.II.estimates,
		as.factor=FALSE, ...){
	
	if(as.factor){
		for(o in outcome.variable){
			rds.data[[o]] <- as.factor(rds.data[[o]])
		}
	}
	
	f <- function(v) cumulative.estimate(rds.data,v,est.func,...)
	ests <- lapply(outcome.variable,f)
	
	make.plot <- function(i){
		#for R CMD check
		Var1 <- Var2 <- value <- NULL
		
		e <- ests[[i]]
		nm <- outcome.variable[i]
		if(ncol(e)==2){
			e1 <- e[,2,drop=FALSE]
			attr(e1,"n") <- attr(e,"n")
			e <- e1
			nm <- paste0(outcome.variable[i],"=",colnames(e)[1])
			rds.data[[outcome.variable[i]]] <- as.factor(rds.data[[outcome.variable[i]]])
		}
		if(ncol(e)>1){
			rownames(e) <- attr(e,"n")
			dat <- melt(e)
			datl <- melt(e[nrow(e),,drop=FALSE])
			p <- ggplot(dat) + 
					geom_line(aes(x=Var1,color=as.factor(Var2),y=value)) + 
					scale_color_hue(nm) +
					ylab("Estimate") + 
					xlab("# of Observations") + 
					scale_y_continuous(limits=c(0,1)) +
					theme_bw()
			p <- p + geom_hline(data=datl,
					aes(yintercept=value,color=as.factor(Var2)),linetype=2,alpha=.5)
			p
		}else{
			dat <- data.frame(value=e[,1],Var1=attr(e,"n"))
			datl <- dat[nrow(dat),,drop=FALSE]
			v <- rds.data[[outcome.variable[i]]]
			rng <- if(!is.numeric(v)) c(0,1) else range(v,na.rm=TRUE)
			p <- ggplot(dat) + 
					geom_line(aes(x=Var1,y=value)) + 
					ylab(paste("Estimated",nm)) + 
					xlab("# of Observations") +
					scale_y_continuous(limits=rng) +
					theme_bw()
			p <- p + geom_hline(data=datl,
					aes(yintercept=value),linetype=2,alpha=.5)
			p
		}
		return(p + ggtitle(paste("Convergence plot of",nm)))
	}
 	plots <- lapply(1:length(outcome.variable),make.plot)
	do.call(grid.arrange,plots)
}

#' Bottleneck Plot
#' @param rds.data An rds.data.frame.
#' @param outcome.variable A character vector of outcome variables.
#' @param est.func A function taking rds.data and outcome.variable as parameters and
#' 	returning an rds.weighted.estimate object.
#' @param ... additional parameters for est.func.
#' @param as.factor Convert all outcome variables to factors
#' @references 
#' Krista J. Gile, Lisa G. Johnston, Matthew J. Salganik \emph{Diagnostics for Respondent-driven Sampling} eprint arXiv:1209.6254, 2012 
#' @examples
#' data(fauxmadrona)
#' bottleneck.plot(fauxmadrona,"disease")
#' @export
bottleneck.plot <- function(rds.data, outcome.variable, est.func=RDS.II.estimates,
		as.factor=FALSE, ...){
	#For R CMD check
	n <- value <- Seed <- NULL
	
	for(o in outcome.variable){
		if(as.factor || is.character(rds.data[[o]])){		
			rds.data[[o]] <- as.factor(rds.data[[o]])
		}
	}
	
	f <- function(v,dat){
		est <- cumulative.estimate(dat,v,est.func,...)
		n <- attr(est,"n")
		if(ncol(est)==1){
			colnames(est) <- v
			rownames(est) <- n
			est <- list(est)
		}else if(ncol(est)==2){
			nl <- colnames(est)[2]
			est <- est[,2,drop=FALSE]
			attr(est,"n") <- n
			colnames(est) <- paste0(v,"=",nl)
			est <- list(est)
		}else{
			est <- lapply(1:ncol(est),function(i){
						e <- est[,i,drop=FALSE]
						nl <- colnames(e)
						attr(e,"n") <- n
						colnames(e) <- paste0(v,"=",nl)
						e
					})
		}
		est
	}
	
	seeds <- get.seed.id(rds.data)
	sids <- unique(seeds)
	ls <- list()
	for(i in 1:length(sids)){
		dat <- rds.data[seeds == sids[i],]
		if(nrow(dat)==0)
			next
		res <- NULL
		nres <- NULL
		for(v in outcome.variable){
			tmp <- f(v,dat)
			for(j in 1:length(tmp)){
				res <- cbind(res,tmp[[j]])
				nres <- cbind(nres,attr(tmp[[j]],"n"))
			}
		}
		res <- melt(res)[-1]
		res$n <- melt(nres)[[3]]
		res$seed <- sids[i]
		ls[[i]] <- res
		names(ls)[i] <- sids[i]
	}
	result <- Reduce(function(a,b){
				if(is.null(a))
					return(b)
				if(is.null(b))
					return(a)
				rbind(a,b)
			},ls,init=NULL)
	result$Seed <- as.factor(result$seed)
	ggplot(result,aes(x=n,y=value,color=Seed)) + 
			geom_line() + 
			facet_wrap(~Var2,scales="free_y") + 
			theme_bw() +
			ylab("Estimate") +
			xlab("# of Observations")
}

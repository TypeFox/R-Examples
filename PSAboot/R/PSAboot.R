#' Returns a vector with the default methods used by \code{\link{PSAboot}}.
#' 
#' The current default methods are:
#' \describe{
#' \item{Stratification}{\code{\link{boot.strata}}}
#' \item{ctree}{\code{\link{boot.ctree}}}
#' \item{rpart}{\code{\link{boot.rpart}}}
#' \item{Matching}{\code{\link{boot.matching}}}
#' \item{MatchIt}{\code{\link{boot.matchit}}}
#' }
#' 
#' The default methods can be changed by setting the \code{PSAboot.methods} option
#' using \code{options('PSAboot.methods'=c(...))} where \code{...} is a named
#' list of functions.
#' 
#' @export
#' @return a vector of methods for use by \code{\link{PSAboot}}
getPSAbootMethods <- function() {
	methods <- getOption('PSAboot.methods', default=c(
		'Stratification'=boot.strata,
		'ctree'=boot.ctree,
		'rpart'=boot.rpart,
		'Matching'=boot.matching,
		'MatchIt'=boot.matchit))
	invisible(methods)
}

#' Bootstraping for propensity score analysis
#' 
#' Bootstrapping has become a popular resampling method for estimating sampling
#' distributions. And propensity score analysis (PSA) has become popular for estimating
#' causal effects in observational studies. This function implements bootstrapping
#' specifically for PSA. Like typical bootstrapping methods, this function estiamtes
#' treatment effects for M random samples. However, unlike typical bootstrap methods,
#' this function allows for separate sample sizes for treatment and control units.
#' That is, under certain circumstances (e.g. when the ratio of treatment-to-control
#' units is large) bootstrapping only the control units may be desirable. Additionally,
#' this function provides a framework to use multiple PSA methods for each bootstrap
#' sample.
#' 
#' @param Tr numeric (0 or 1) or logical vector of treatment indicators. 
#' @param Y vector of outcome varaible.
#' @param X matrix or data frame of covariates used to estimate the propensity scores.
#' @param M number of bootstrap samples to generate.
#' @param formu formula used for estimating propensity scores. The default is to use
#'        all covariates in \code{X}.
#' @param control.ratio the ratio of control units to sample relative to the treatment units.
#' @param control.sample.size the size of each bootstrap sample of control units.
#' @param control.replace whether to use replacement when sampling from control units.
#' @param treated.sample.size the size of each bootstrap sample of treatment units. The
#'        default uses all treatment units for each boostrap sample.
#' @param treated.replace whether to use replacement when sampling from treated units.
#' @param methods a named vector of functions for each PSA method to use.
#' @param seed random seed. Each iteration, i, will use a seed of \code{seed + i}.
#' @param parallel whether to run the bootstrap samples in parallel.
#' @param ... other parameters passed to \code{\link{Match}} and \code{\link{psa.strata}}
#' @return a list with following elements:
#' 		  \describe{
#' 		  \item{overall.summary}{Data frame with the results using the complete
#' 		                         dateset (i.e. unboostrapped results).}
#' 		  \item{overall.details}{Objects returned from each method for complete dataset.}
#' 		  \item{pooled.summary}{Data frame with results of each boostrap sample.}
#' 		  \item{pooled.details}{List of objects returned from each method for each 
#' 		        boostrap sample.}
#' 		  \item{control.sample.size}{sample size used for control units.}
#' 		  \item{treated.sample.size}{sample size used for treated units.}
#' 		  \item{control.replace}{whether control units were sampled with replacement.}
#' 		  \item{treated.replace}{whether treated units were sampled with replacement.}
#' 		  \item{Tr}{vector of treatment assignment.}
#' 		  \item{Y}{vector out outcome.}
#' 		  \item{X}{matrix or data frame of covariates.}
#' 		  \item{M}{number of bootstrap samples.}
#' 		  }
#' @seealso getPSAbootMethods
#' @export
PSAboot <- function(Tr, Y, X, M=100,
					formu=as.formula(paste0('treat ~ ', paste0(names(X), collapse=' + '))),
					control.ratio=5,
					control.sample.size=min(control.ratio*min(table(Tr)),
											max(table(Tr))),
					control.replace=TRUE,
					treated.sample.size=min(table(Tr)),
					treated.replace=TRUE,
					methods=getPSAbootMethods(),
					parallel=TRUE,
					seed=NULL,
					...) {
	if(length(Tr) != length(Y) & length(Y) != nrow(X)) {
		stop('The length of Tr and Y and the number of rows of X must be the same!')
	}
	if('factor' %in% class(Tr)) {
		groups <- levels(Tr)
		message(paste0('Using ', groups[2], ' as treatment group.'))
	} else {
		groups <- 0:1
	}
	index.control <- which(Tr == groups[1])
	index.treated <- which(Tr == groups[2])
	if(!control.replace & control.sample.size > length(index.control)) { 
		stop('Sample size cannot be larger than the number of control units. 
			 Try a smaller control.ratio or specify control.replace=TRUE.')
	}
	
	message(paste0(
		   M, ' bootstrap samples using ', length(methods), ' methods.\n',
		   'Bootstrap sample sizes:\n',
		   '   Treated=', length(index.treated), ' (', 
		   round(treated.sample.size / length(index.treated) * 100), '%) with', 
		   ifelse(treated.replace, '', 'out'), ' replacement.\n',
		   '   Control=', length(index.control), ' (', 
		   round(control.sample.size / length(index.control) * 100), '%) with', 
		   ifelse(control.replace, '', 'out'), ' replacement.'
	))
	
	tmp.formu <- update.formula(formu, 'treat ~ .')
	vars <- all.vars(tmp.formu)
	vars <- vars[2:length(vars)] #Keep only the right hand variables
	X <- X[,vars] #This will ensure that X contains only the covariates
	if('factor' %in% sapply(X, class)) {
		X.trans <- cv.trans.psa(X)[[1]]
	} else {
		X.trans <- X
	}
	complete.summary <- data.frame()
	complete.details <- list()
	for(m in seq_along(methods)) {
		n <- names(methods)[[m]]
		f <- methods[[m]]
		r <- f(Tr=Tr, Y=Y, X=X, X.trans=X.trans, formu=formu, ...)
		complete.details[[paste0('summary.', n)]] <- r$summary
		complete.details[[paste0('details.', n)]] <- r$details
		complete.details[[paste0('balance.', n)]] <- r$balance
		complete.summary <- rbind(complete.summary, data.frame(
			method=n, 
			estimate=unname(r$summary['estimate']),
			ci.min=unname(r$summary['ci.min']),
			ci.max=unname(r$summary['ci.max']),
			stringsAsFactors=FALSE))
	}
	
	pb <- txtProgressBar(1,M,style=3)
	
	bootfun <- function(i) {
		if(!is.null(seed)) { set.seed(seed + i) }
		index.control.sample <- sample(index.control, size=control.sample.size, 
									   replace=control.replace)
		index.treated.sample <- sample(index.treated, size=treated.sample.size, 
									   replace=treated.replace)
		rows <- c(index.treated.sample, index.control.sample)
		result <- list()
		result$index.control <- index.control.sample
		result$index.treated <- index.treated.sample
		result[['summary']] <- data.frame()
		for(m in seq_along(methods)) {
			n <- names(methods)[[m]]
			f <- methods[[m]]
			tryCatch({
				r <- f(Tr=Tr[rows], Y=Y[rows], X=X[rows,], X.trans=X.trans[rows,], 
					   formu=formu, ...)
				result[[paste0('summary.', n)]] <- r$summary
				result[[paste0('details.', n)]] <- r$details
				result[[paste0('balance.', n)]] <- r$balance
				result[['summary']] <- rbind(result[['summary']], data.frame(
					Method=n, 
					estimate=unname(r$summary['estimate']),
					ci.min=unname(r$summary['ci.min']),
					ci.max=unname(r$summary['ci.max']),
					stringsAsFactors=FALSE))
			}, error=function(e) { 
				warning(paste0('Error occurred during iteration ', i,
							   ' for ', n, ' method: ', e))
			})
		}
		setTxtProgressBar(pb, i); flush.console()
		return(result)
	}
	
	if(parallel) {
		tmp <- mclapply(seq_len(M), FUN=bootfun)
	} else {
		tmp <- lapply(seq_len(M), FUN=bootfun)
	}
	
	close(pb)
	
	summary <- data.frame(iter=rep(1:M, each=length(methods)),
						  method=rep(names(methods), M),
						  estimate=rep(as.numeric(NA), M * length(methods)),
						  ci.min=rep(as.numeric(NA), M * length(methods)),
						  ci.max=rep(as.numeric(NA), M * length(methods)),
						  stringsAsFactors=FALSE)
	cols <- c('estimate','ci.min','ci.max')
	for(i in seq_along(tmp)) {
		sum <- tmp[[i]]$summary
		for(j in 1:nrow(sum)) {
			summary[summary$iter == i & summary$method==sum[j,]$Method, 
					cols] <- as.numeric((sum[j,cols]))
		}
	}
	r <- list(pooled.summary=summary,
			  pooled.details=tmp,
			  complete.summary=complete.summary,
			  complete.details=complete.details, 
			  Y=Y, Tr=Tr, X=X, M=M, seed=seed,
			  control.sample.size=control.sample.size,
			  treated.sample.size=treated.sample.size,
			  control.replace=control.replace,
			  treated.replace=treated.replace)
	class(r) <- "PSAboot"
	return(r)
}

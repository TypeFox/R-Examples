#' Estimate covariate effect sizes before and after propensity score adjustment.
#' 
#' @param covariates frame or matrix of covariates.
#' @param treatment vector of treatment indicators.
#' @param level2 vector indicating level 2 membership.
#' @param strata strata indicators.
#' @param abs if TRUE absolute values of effect sizes will be plotted.
#' @export
covariate.balance <- function(covariates, treatment, level2, strata, abs=TRUE) {
	#Recode factors. First we'll covert logicals and factors with two levels to integers
	for(i in 1:ncol(covariates)) {
		if(class(covariates[,i]) == 'logical') {
			covariates[,i] <- as.integer(covariates[,i])
		} else if(class(covariates[,i]) == 'factor' & length(levels(covariates[,i])) == 2) {
			covariates[,i] <- as.integer(covariates[,i])
		}
	}
	if('factor' %in% sapply(covariates, class)) {
		#Convert remaining factors using cv.trans.psa from PSAgraphics
		#covariates <- cv.trans.psa(covariates)[[1]]
		tmp.covariates <- covariates[,sapply(covariates, class) == 'factor', drop=FALSE]
		covariates <- covariates[,sapply(covariates, class) != 'factor', drop=FALSE]
		covariates <- as.data.frame(rbind(model.matrix(~ ., tmp.covariates)))
	}
	
	results <- data.frame(row.names=names(covariates), 
						  es.adj=rep(as.numeric(NA), ncol(covariates)), 
						  es.adj.wtd=rep(as.numeric(NA), ncol(covariates)),
						  es.unadj=rep(as.numeric(NA), ncol(covariates)),
						  stringsAsFactors=FALSE)
	strata.es <- list()
	for(i in names(covariates)) {
		rows <- !is.na(covariates[,i])
		adj <- data.frame()
		for(l in unique(level2)) {
			lvl.rows <- rows & level2 == l
			tab <- cbind(
				as.data.frame(tapply(covariates[lvl.rows,i], 
									 list(strata[lvl.rows], treatment[lvl.rows]), mean)),
				as.data.frame(tapply(covariates[lvl.rows,i], 
									 list(strata[lvl.rows]), length)),
				as.data.frame(tapply(covariates[lvl.rows,i],
									 list(strata[lvl.rows]), sd))
			)
			names(tab)[3:4] <- c('n','sd')
			tab$level2 <- l
			tab$strata <- row.names(tab)
			adj <- rbind(adj, tab)
		}
		adj <- na.omit(adj)
		adj$diff <- adj[,1] - adj[,2]
		adj$es <- adj$diff / adj$sd
		nans <- is.nan(adj$es)
		if(any(nans)) {
			adj[is.nan(adj$es),]$es <- 0 # If there is perfect balance then the st dev will be 0
		}
		es.wtd <- adj$es * adj$n
		results[i,]$es.adj <- mean(adj$es)
		results[i,]$es.adj.wtd <- mean(es.wtd) / sum(adj$n)
		results[i,]$es.unadj <- ( mean(covariates[treatment == names(adj)[1],i], na.rm=TRUE) - 
								  mean(covariates[treatment == names(adj)[2],i], na.rm=TRUE) ) / 
			                      sd(covariates[,i], na.rm=TRUE)
		strata.es[[i]] <- adj
	}
	if(abs) {
		results$es.adj <- abs(results$es.adj)
		results$es.unadj <- abs(results$es.unadj)
	}
	results <- cbind(covariate=row.names(results), results)
	row.names(results) <- 1:nrow(results)
	ret <- list(effects=results, strata.effects=strata.es)
	class(ret) <- 'covariate.balance'
	return(ret)
}

#' Returns the overall effects as a data frame.
#' 
#' @param x results of \code{\link{covariate.balance}}.
#' @param row.names unused.
#' @param optional unused.
#' @param ... unused
#' @return a data frame with overall covariate effects before and after adjustment.
#' @method as.data.frame covariate.balance
#' @export
as.data.frame.covariate.balance <- function(x, row.names=NULL, optional=FALSE, ...) {
	return(x$effects)
}

#' Prints the overall effects before and after propensity score adjustment.
#' 
#' @param x results of \code{\link{covariate.balance}}.
#' @param ... unused.
#' @method print covariate.balance
#' @export
print.covariate.balance <- function(x, ...) {
	print(x$effects)
}

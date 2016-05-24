#' Forward search algorithm for outlier detection
#'
#' The forward search algorithm begins by selecting a homogeneous subset
#' of cases based on a maximum likelihood criteria and continues to add individual
#' cases at each iteration given an acceptance criteria. By default the function
#' will add cases that contribute most to the likelihood function and that have
#' the closest robust Mahalanobis distance, however model implied residuals
#' may be included as well.
#'
#' Note that \code{forward.search} is not limited to confirmatory factor analysis and
#' can apply to nearly any model being studied
#' where detection of influential observations is important.
#'
#' @aliases forward.search
#' @param data matrix or data.frame
#' @param model if a single numeric number declares number of factors to extract in
#'   exploratory factor analysis. If \code{class(model)} is a sem (semmod), or lavaan (character),
#'   then a confirmatory approach is performed instead
#' @param criteria character strings indicating the forward search method
#'   Can contain \code{'GOF'} for goodness of fit distance, \code{'mah'} for Mahalanobis
#'   distance, or \code{'res'} for model implied residuals
#' @param n.subsets a scalar indicating how many samples to draw to find
#'   a homogeneous starting base group
#' @param p.base proportion of sample size to use as the base group
#' @param print.messages logical; print how many iterations are remaining?
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#'   \code{\link{gCD}}, \code{\link{GOF}}, \code{\link{LD}},
#'   \code{\link{robustMD}}, \code{\link{setCluster}}
#' @keywords forward.search
#' @export forward.search
#' @examples
#'
#' \dontrun{
#'
#' #run all internal gCD and GOF functions using multiple cores
#' setCluster()
#'
#' #Exploratory
#' nfact <- 3
#' (FS <- forward.search(holzinger, nfact))
#' (FS.outlier <- forward.search(holzinger.outlier, nfact))
#' plot(FS)
#' plot(FS.outlier)
#'
#' #Confirmatory with sem
#' model <- sem::specifyModel()
#'	  F1 -> Remndrs,    lam11
#' 	  F1 -> SntComp,    lam21
#' 	  F1 -> WrdMean,    lam31
#' 	  F2 -> MissNum,    lam41
#' 	  F2 -> MxdArit,    lam52
#' 	  F2 -> OddWrds,    lam62
#' 	  F3 -> Boots,      lam73
#'	  F3 -> Gloves,     lam83
#' 	  F3 -> Hatchts,    lam93
#' 	  F1 <-> F1,   NA,     1
#' 	  F2 <-> F2,   NA,     1
#' 	  F3 <-> F3,   NA,     1
#'
#'
#' (FS <- forward.search(holzinger, model))
#' (FS.outlier <- forward.search(holzinger.outlier, model))
#' plot(FS)
#' plot(FS.outlier)
#'
#' #Confirmatory with lavaan
#' model <- 'F1 =~  Remndrs + SntComp + WrdMean
#' F2 =~ MissNum + MxdArit + OddWrds
#' F3 =~ Boots + Gloves + Hatchts'
#'
#' (FS <- forward.search(holzinger, model))
#' (FS.outlier <- forward.search(holzinger.outlier, model))
#' plot(FS)
#' plot(FS.outlier)
#'
#'
#' }
forward.search <- function(data, model, criteria = c('GOF', 'mah'),
	n.subsets = 1000, p.base= .4, print.messages = TRUE, ...)
{
    if(any(is.na(data)))
        stop('All routines require complete datasets (no NA\'s) so that the search
             gives meaninful results.')
	N <- nrow(data)
	p <- ncol(data)
	ID <- 1:N
	Samples <- matrix(0, floor(p.base*N), n.subsets)
	for(i in 1:n.subsets)
		Samples[ ,i] <- sample(1:N, floor(p.base*N))
	if(is.numeric(model)){
		STATISTICS <- myApply(matrix(1:n.subsets), 1, function(i, data, Samples, model){
            ret <- try(factanal(data[Samples[ ,i], ], model)$STATISTIC, TRUE)
            if(is(ret, 'try-error')) ret <- Inf
		    return(ret)
		}, data=data, Samples=Samples, model=model)
		orgbaseID <- baseID <- Samples[ ,(min(STATISTICS) == STATISTICS)]
		nbaseID <- setdiff(ID, baseID)
		basedata <- data[baseID, ]
		basemodels <- list()
		orderentered <- c()
		for (LOOP in 1:length(nbaseID)){
			tmpcor <- cor(basedata)
			basemodels[[LOOP]] <- mlfact(tmpcor, model)
			basemodels[[LOOP]]$N <- nrow(basedata)
			basemodels[[LOOP]]$R <- tmpcor
 			stat <- c()
			RANK <- rep(0, length(nbaseID))
			if(any(criteria == 'GOF')){
			    stat <- myApply(matrix(1:length(nbaseID)), 1, function(j, basedata, nbaseID, model){
			        ret <- try(mlfact(cor(rbind(basedata, data[nbaseID[j], ])), model)$value, TRUE)
			        if(is(ret, 'try-error')) ret <- Inf
			        return(ret)
			    }, basedata=basedata, nbaseID=nbaseID, model=model)
				RANK <- RANK + rank(stat)
			}
			if(any(criteria == 'mah')){
				stat <- mahalanobis(data[nbaseID, ], colMeans(data[baseID, ]),
					cov(data[baseID, ]))
				RANK <- RANK + rank(stat)
			}
			if(any(criteria == 'res')){
				stat <- c()
				for(j in 1:length(nbaseID)){
					tmp <- obs.resid(rbind(basedata, data[nbaseID[j], ]), model)
					stat[j] <- tmp$std_res[nrow(basedata)+1, ] %*%
						tmp$std_res[nrow(basedata)+1, ]
				}
				RANK <- RANK + rank(stat)
			}
			RANK <- rank(RANK)
			newID <- nbaseID[min(RANK) == RANK]
			if(length(newID) > 1) newID <- newID[1]
			orderentered <- c(orderentered, newID)
			baseID <- c(baseID, newID)
			nbaseID <- setdiff(ID, baseID)
			basedata <- data[baseID, ]
			if(print.messages) {
				cat('Remaining iterations:', length(nbaseID), '\n')
				flush.console()
			}
		}
		basemodels[[LOOP+1]] <- mlfact(cor(data), model)
		basemodels[[LOOP+1]]$N <- nrow(data)
		basemodels[[LOOP+1]]$R <- cor(data)
		GOFstat <- RMR <- Cooksstat <- c()
		for(i in 1:(length(basemodels)-1)){
			GOFstat[i] <- basemodels[[i]]$value * (length(orgbaseID) + i - 1)
			theta <- basemodels[[i]]$par
			hess <- basemodels[[i]]$hessian
			theta2 <- basemodels[[i+1]]$par
			Cooksstat[i] <- (theta-theta2) %*% solve(hess) %*% (theta-theta2)
			Rhat <- basemodels[[i]]$loadings %*% t(basemodels[[i]]$loadings)
			diag(Rhat) <- 1
			RMR[i] <- sqrt(2*sum(((basemodels[[i]]$R - Rhat)^2) /
				(ncol(Rhat)*(ncol(Rhat) + 1))))
		}
		Cooksstat <- c(NA, Cooksstat)
		orderentered <- c(NA, orderentered)
		GOFstat[i+1] <- basemodels[[i+1]]$value * N
		Rhat <- basemodels[[i+1]]$loadings %*% t(basemodels[[i+1]]$loadings)
		diag(Rhat) <- 1
		RMR[i+1] <- sqrt(2*sum(((basemodels[[i+1]]$R - Rhat)^2) /
			(ncol(Rhat)*(ncol(Rhat) + 1))))
		ret <- list(GOF=GOFstat, RMR=RMR, gCD=Cooksstat, ord=orderentered)
	} else if(class(model) == "semmod"){
	    sampleCov <- cov(data)
	    STATISTICS <- myApply(matrix(1:n.subsets), 1, function(i, data, Samples, model){
            tmpdat <- data[Samples[ ,i], ]
	        samplesemMod <- try(sem::sem(model, cov(tmpdat), nrow(tmpdat), ...), TRUE)
            if(is(samplesemMod, 'try-error')) return(Inf)
            ret <- samplesemMod$criterion * (samplesemMod$N - 1)
            return(ifelse(ret < 0, Inf, ret))
	    }, data=data, Samples=Samples, model=model)
        orgbaseID <- baseID <- Samples[ ,(min(STATISTICS) == STATISTICS)]
	    nbaseID <- setdiff(ID, baseID)
	    basedata <- data[baseID, ]
	    stat <- c()
	    basemodels <- list()
	    orderentered <- c()
	    for (LOOP in 1:length(nbaseID)){
	        tmpcov <- cov(basedata)
	        basemodels[[LOOP]] <- sem::sem(model, tmpcov, nrow(basedata), ...)
	        RANK <- rep(0, length(nbaseID))
			if(any(criteria == 'GOF')){
                stat <- myApply(matrix(1:length(nbaseID)), 1, function(j, basedata, nbaseID, model){
                    tmpcov <- cov(rbind(basedata, data[nbaseID[j], ]))
                    tmpmod <- try(sem::sem(model, tmpcov, nrow(basedata) + 1, ...), TRUE)
                    if(is(tmpmod, 'try-error')) return(Inf)
                    ret <- tmpmod$criterion * (tmpmod$N - 1)
                    return(ifelse(ret < 0, Inf, ret))
                }, basedata=basedata, nbaseID=nbaseID, model=model)
				RANK <- RANK + rank(stat)
			}
			if(any(criteria == 'mah')){
				stat <- mahalanobis(data[nbaseID, ], colMeans(data[baseID, ]),
					cov(data[baseID, ]))
				RANK <- RANK + rank(stat)
			}
			if(any(criteria == 'res')){
				stat <- c()
				for(j in 1:length(nbaseID)){
					tmp <- obs.resid(rbind(basedata, data[nbaseID[j], ]), model)
					stat[j] <- tmp$std_res[nrow(basedata)+1, ] %*%
						tmp$std_res[nrow(basedata)+1, ]
				}
				RANK <- RANK + rank(stat)
			}
	    	RANK <- rank(RANK)
	    	newID <- nbaseID[min(RANK) == RANK]
	    	if(length(newID) > 1) newID <- newID[1]
	    	orderentered <- c(orderentered, newID)
	    	baseID <- c(baseID, newID)
	    	nbaseID <- setdiff(ID, baseID)
	    	basedata <- data[baseID, ]
	    	if(print.messages){
	    		cat('Remaining iterations:', length(nbaseID), '\n')
	    		flush.console()
	    	}
	    }
	    tmpcov <- cov(data)
	    basemodels[[LOOP+1]] <- sem::sem(model, tmpcov, N, ...)
	    GOFstat <- RMR <- Cooksstat <- c()
	    for(i in 1:(length(basemodels)-1)){
	    	GOFstat[i] <- basemodels[[i]]$criterion * (basemodels[[i]]$N - 1)
	    	theta <- basemodels[[i]]$coeff
	    	vcov <- basemodels[[i]]$vcov
	    	theta2 <- basemodels[[i+1]]$coeff
	    	Cooksstat[i] <- (theta-theta2) %*% vcov %*% (theta-theta2)
	    	Chat <- basemodels[[i]]$C
	    	C <- basemodels[[i]]$S
	    	RMR[i] <- sqrt(2*sum(((C - Chat)^2) /
	    		(ncol(C)*(ncol(C) + 1))))
	    }
	    Cooksstat <- c(NA, Cooksstat)
	    orderentered <- c(NA, orderentered)
	    GOFstat[i+1] <- basemodels[[i+1]]$criterion * (basemodels[[i+1]]$N - 1)
	    Chat <- basemodels[[i+1]]$C
	    C <- basemodels[[i+1]]$S
	    RMR[i+1] <- sqrt(2*sum(((C - Chat)^2) /	(ncol(C)*(ncol(C) + 1))))
	    ret <- list(GOF=GOFstat, RMR=RMR, gCD=Cooksstat, ord=orderentered)
	} else if(class(model) == "character"){
	    STATISTICS <- myApply(matrix(1:n.subsets), 1, function(i, data, Samples, model){
	        tmpdat <- data[Samples[ ,i], ]
	        samplesemMod <- try(lavaan::sem(model, data=tmpdat, ...), TRUE)
            if(is(samplesemMod, 'try-error')) return(Inf)
            ret <- samplesemMod@Fit@test[[1]]$stat
	        return(ifelse(is.na(ret), Inf, ret))
	    }, data=data, Samples=Samples, model=model)
	    orgbaseID <- baseID <- Samples[ ,(min(STATISTICS) == STATISTICS)]
	    nbaseID <- setdiff(ID, baseID)
	    basedata <- data[baseID, ]
	    basemodels <- list()
	    orderentered <- c()
	    for (LOOP in 1:length(nbaseID)){
	        basemodels[[LOOP]] <- lavaan::sem(model, data=basedata, ...)
	        stat <- c()
	        RANK <- rep(0, length(nbaseID))
	        if(any(criteria == 'GOF')){
	            stat <- myApply(matrix(1:length(nbaseID)), 1, function(j, basedata, nbaseID, model){
	                tmpdat <- rbind(basedata, data[nbaseID[j], ])
	                tmpmod <- try(lavaan::sem(model, data=tmpdat, ...), TRUE)
	                if(is(tmpmod, 'try-error')) return(Inf)
	                ret <- tmpmod@Fit@test[[1]]$stat
	                return(ifelse(is.na(ret), Inf, ret))
	            }, basedata=basedata, nbaseID=nbaseID, model=model)
	            RANK <- RANK + rank(stat)
	        }
	        if(any(criteria == 'mah')){
	            stat <- mahalanobis(data[nbaseID, ], colMeans(data[baseID, ]),
	                                cov(data[baseID, ]))
	            RANK <- RANK + rank(stat)
	        }
	        if(any(criteria == 'res')){
	            stat <- c()
	            for(j in 1:length(nbaseID)){
	                tmp <- obs.resid(rbind(basedata, data[nbaseID[j], ]), model)
	                stat[j] <- tmp$std_res[nrow(basedata)+1, ] %*%
	                    tmp$std_res[nrow(basedata)+1, ]
	            }
	            RANK <- RANK + rank(stat)
	        }
	        RANK <- rank(RANK)
	        newID <- nbaseID[min(RANK) == RANK]
	        if(length(newID) > 1) newID <- newID[1]
	        orderentered <- c(orderentered, newID)
	        baseID <- c(baseID, newID)
	        nbaseID <- setdiff(ID, baseID)
	        basedata <- data[baseID, ]
	        if(print.messages){
	            cat('Remaining iterations:', length(nbaseID), '\n')
	            flush.console()
	        }
	    }
	    basemodels[[LOOP+1]] <- lavaan::sem(model, data=data,  ...)
	    GOFstat <- RMR <- Cooksstat <- c()
	    for(i in 1:(length(basemodels)-1)){
	        GOFstat[i] <- basemodels[[i]]@Fit@test[[1]]$stat
	        theta <- basemodels[[i]]@Fit@x
	        vcov <- lavaan::vcov(basemodels[[i]])
	        theta2 <- basemodels[[i+1]]@Fit@x
	        Cooksstat[i] <- (theta-theta2) %*% vcov %*% (theta-theta2)
	        Chat <- basemodels[[i]]@Fit@Sigma.hat[[1]]
	        C <- basemodels[[i]]@SampleStats@cov[[1]]
	        RMR[i] <- sqrt(2*sum(((C - Chat)^2) /
	                                 (ncol(C)*(ncol(C) + 1))))
	    }
	    Cooksstat <- c(NA, Cooksstat)
	    orderentered <- c(NA, orderentered)
	    GOFstat[i+1] <- basemodels[[i+1]]@Fit@test[[1]]$stat
	    Chat <- basemodels[[i+1]]@Fit@Sigma.hat[[1]]
	    C <- basemodels[[i+1]]@SampleStats@cov[[1]]
	    RMR[i+1] <- sqrt(2*sum(((C - Chat)^2) /	(ncol(C)*(ncol(C) + 1))))
	    ret <- list(GOF=GOFstat, RMR=RMR, gCD=Cooksstat, ord=orderentered)
	} else {
	    stop('model class not supported')
	}

	class(ret) <- 'forward.search'
	ret
}

#' @rdname forward.search
#' @param x an object of class \code{forward.search}
#' @param ncases number of final cases to print in the sequence
#' @param stat type of statistic to use. Could be 'GOF', 'RMR', or 'gCD' for
#'   the model chi squared value, root mean square residual, or generalized Cook's distance,
#'   respectively
#' @param ... additional parameters to be passed
#' @export
print.forward.search <- function(x, ncases = 10, stat = 'GOF', ...)
{
	if(stat == 'GOF') ret <- x$GOF
	if(stat == 'RMR') ret <- x$RMR
	if(stat == 'gCD') ret <- x$gCD
	names(ret) <- x$ord
	print(tail(ret, ncases))
	return(invisible())
}

#' @rdname forward.search
#' @param y a \code{null} value ignored by \code{plot}
#' @param main the main title of the plot
#' @param type type of plot to use, default displays points and lines
#' @param ylab the y label of the plot
#' @export
plot.forward.search <- function(x, y = NULL, stat = 'GOF', main = 'Forward Search',
	type = c('p','h'), ylab = 'obs.resid', ...)
{
	id <- x$ord
	Input <- 1:length(id)
	if(stat == 'GOF') stat2 <- x$GOF
	if(stat == 'RMR') stat2 <- x$RMR
	if(stat == 'gCD') stat2 <- x$gCD
	dat <- data.frame(stat2,Input,id)
	ret <- lattice::xyplot(stat2~Input, dat, type=type, main=main, ylab=stat,
		xlab='Input Step', ...)
	return(ret)
}


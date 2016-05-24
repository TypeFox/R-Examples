bic.mlogit <- function (f, data, choices=NULL, base.choice=1,
							varying = NULL, sep='.', approx=TRUE, 
							include.intercepts=TRUE, 
							verbose=FALSE, ...
						) {
	cl <- match.call()
		
	bin.logit <- mlogit2logit(f, data=data, choices=choices, base.choice=base.choice, 
								varying=varying, sep=sep)
	
	mnl <- list()
	mnl[['specification']] <- mnl.spec(f, data=data, choices=choices, base.choice=base.choice, 
										varying=varying, sep=sep)
	mnl[['zcols']] <- bin.logit$zcols
	mnl[['choice.main.intercept']] <- bin.logit$choice.main.intercept
				
	which.y <- which.max(mnl$specification$response == colnames(data))
	mnl[['x']] <- data[-which.y]
	mnl[['y']] <- data[,which.y]
	varying.tmp <- rep(FALSE, ncol(data))
	varying.tmp[varying] <- TRUE
	varying.tmp <- varying.tmp[-which.y]
	mnl[['varying']] <- cumsum(varying.tmp)

	bma.res <- bic.glm.bg(f=bin.logit$formula, data=bin.logit$data, mnl=mnl, 
									glm.family='binomial', approx=approx, 
									include.intercepts=include.intercepts,
									verbose=verbose, ...
							)
	bma.res$bic.glm$call <- cl
	invisible(structure(list(bic.glm=bma.res$bic.glm, bin.logit=bin.logit, spec=mnl[['specification']],
							bma.specifications=bma.res$specifications, approx=approx),
				class='bic.mlogit'))
}

summary.bic.mlogit <- function(object, ...) {
	summary(object$bic.glm, ...)
	cat('\nMNL specification:\n')
	cat('==================')
	summary(object$spec, ...)
}

imageplot.mlogit <- function(x, ...) imageplot.bma(x$bic.glm, ...)

plot.bic.mlogit <- function(x, ...) plot(x$bic.glm, ...)

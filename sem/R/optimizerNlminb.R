# last modified 2011-07-30

optimizerNlminb <- function(start, objective=objectiveML, 
		gradient=TRUE, maxiter, debug, par.size, model.description, warn, ...){
	with(model.description, {
				obj <- objective(gradient=gradient)
				objective <- obj$objective
				grad <- if (gradient) obj$gradient else NULL
				if (!warn) save.warn <- options(warn=-1)
				res <- nlminb(start, objective, gradient=grad, model.description=model.description, 
						control=list(trace=if(debug) 1 else 0, iter.max=maxiter, ...))
				if (!warn) options(save.warn)
				result <- list()
				result$convergence <- res$convergence == 0
				result$iterations <- res$iterations
				par <- res$par
				names(par) <- param.names
				result$par <- par
				if (!result$convergence)
					warning(paste('Optimization may not have converged; nlminb return code = ',
									res$convergence, '. Consult ?nlminb.\n', sep=""))
				result$criterion <- res$objective
				obj <- objective(par, model.description)
				C <- attr(obj, "C")
				rownames(C) <- colnames(C) <- var.names[observed]
				result$C <- C
				A <- attr(obj, "A")
				rownames(A) <- colnames(A) <- var.names
				result$A <- A
				P <- attr(obj, "P")
				rownames(P) <- colnames(P) <- var.names
				result$P <- P
				class(result) <- "semResult"
				result
			}
	)
}

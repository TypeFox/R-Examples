# last modified 2011-08-10 by J. Fox

optimizerNlm <- function(start, objective=objectiveML, 
	gradient=TRUE, maxiter, debug, par.size, model.description, warn, ...){
	with(model.description, {
			obj <- objective(gradient=gradient)$objective
			typsize <- if (par.size == 'startvalues') abs(start) else rep(1, t)
			if (!warn) save.warn <- options(warn=-1)
			res <- nlm(obj, start, iterlim=maxiter, print.level=if(debug) 2 else 0,
				typsize=typsize, hessian=TRUE, model.description, ...)
			if (!warn) options(save.warn)
			result <- list()
			result$convergence <- res$code <= 2
			result$iterations <- res$iterations
			par <- res$estimate
			names(par) <- param.names
			result$par <- par
			if (!result$convergence)
				warning(paste('Optimization may not have converged; nlm return code = ',
						res$code, '. Consult ?nlm.\n', sep=""))
			vcov <- matrix(NA, t, t)
			qr.hess <- try(qr(res$hessian), silent=TRUE)
			if (class(qr.hess) == "try-error"){
				warning("Could not compute QR decomposition of Hessian.\nOptimization probably did not converge.\n")
			}
			else if (qr.hess$rank < t){
				warning(' singular Hessian: model is probably underidentified.\n')
				which.aliased <- qr.hess$pivot[-(1:qr.hess$rank)]
				result$aliased <- param.names[which.aliased]
			}
			else {
				vcov <- (2/(N - (!raw))) * solve(res$hessian)
				if (any(diag(vcov) < 0)) {
					result$aliased <- param.names[diag(vcov) < 0]
					warning("Negative parameter variances.\nModel may be underidentified.\n")
				}
			}
			colnames(vcov) <- rownames(vcov) <- param.names
			result$vcov <- vcov
			result$criterion <- res$minimum # c(result$obj) - n - log(det(S))
			obj <- obj(par, model.description)
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
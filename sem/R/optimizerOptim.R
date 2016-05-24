# last modified 2011-08-10 by J. Fox

optimizerOptim <- function(start, objective=objectiveML, 
	gradient=TRUE, maxiter, debug, par.size, model.description, warn, method="CG", ...){
	with(model.description, {
			obj <- objective(gradient=gradient)
			grad <- if (gradient) obj$gradient else NULL
			obj <- obj$objective
			typsize <- if (par.size == 'startvalues') abs(start) else rep(1, t)
			if (!warn) save.warn <- options(warn=-1)
			res <- optim(start, obj, gr=grad, hessian=TRUE, control=list(maxit=maxiter, trace=if(debug) 6 else 0,
				parscale=typsize), method=method, ..., model.description)
			if (!warn) options(save.warn)
			result <- list()
			result$convergence <- res$convergence == 0
			result$iterations <- if (res$convergence == 1) maxiter else NA
			par <- res$par
			names(par) <- param.names
			result$par <- par
			if (!result$convergence)
				warning(paste('Optimization may not have converged; optim return code = ',
						res$convergence, '. Consult ?optim.\n', sep=""))
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
			result$criterion <- res$value 
			if (!raw) {
				CC <- diag(diag(S))
				result$chisqNull <- (N - 1) * 
					(sum(diag(S %*% solve(CC))) + log(det(CC)) - log(det(S)) - n)
			}
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
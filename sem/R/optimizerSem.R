# last modified 2012-01-07 by J. Fox
# Modified for Compiled Objective and nlm in C/C++ by Zhenghua Nie.

optimizerSem <- function(start, objective=objectiveML,  
	gradient=TRUE, maxiter, debug, par.size, model.description, warn, ...){
	with(model.description, {
			obj <- objective(gradient=gradient)$objective
			typsize <- if (par.size == 'startvalues') abs(start) else rep(1, t) #we move typesize in csem.cpp.

#			objectiveCompiled <- "objectiveML"
#			if(identical(objective, objectiveCompiledML) || identical(objective, objectiveML)) 
#					objectiveCompiled <- "objectiveML"
#			if(identical(objective, objectiveCompiledGLS) || identical(objective, objectiveGLS))
#					objectiveCompiled <- "objectiveGLS"
				
#			objectiveCompiled <- deparse(substitute(objective))
#			if (!objectiveCompiled %in% c("objectiveML", "objectiveGLS")) stop("optimizerSem requires the objectiveML or objectiveGLS objective function")
#			objectiveCompiled <- "objectiveML"
				
			if(identical(objective, objectiveML)) objectiveCompiled <- "objectiveML"
			else if (identical(objective, objectiveGLS)) objectiveCompiled <- "objectiveGLS"
			else if (identical(objective, objectiveFIML)) objectiveCompiled <- "objectiveFIML"
			else stop("optimizerSem requires the objectiveML or objectiveGLS or objectiveFIML objective function")
			
			if (!warn) save.warn <- options(warn=-1)

			res <- CompiledSolve(model.description=model.description, start=start, objective=objectiveCompiled, gradient=gradient, typsize=typsize, debug=debug, maxiter=maxiter)

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
#			if(identical(objective, objectiveML) || identical(objective, objectiveGLS)){
					C <- res$C
					A <- res$A
					P <- res$P
#			} else {
#					obj <- obj(par, model.description)  
#					C <- attr(obj, "C")
#					A <- attr(obj, "A")
#					P <- attr(obj, "P")
#			}

			rownames(C) <- colnames(C) <- var.names[observed]
			result$C <- C
			rownames(A) <- colnames(A) <- var.names
			result$A <- A
			rownames(P) <- colnames(P) <- var.names
			result$P <- P
			class(result) <- "semResult"
			result
	}
)
}

# last modified 2012-01-06 by J. Fox
# Modified for Compiled code in C/C++ by Zhenghua Nie.

objectiveFIML <- function (gradient = TRUE, hessian=FALSE) 
{
    result <- list(objective = function(par, model.description) {
        with(model.description, {
						 
						 res <- CompiledObjective(par=par, model.description=model.description, objective="objectiveFIML", gradient=gradient, hessian=hessian) 
				
						 f <- res$f
				     C <- res$C
				     A <- res$A
				     P <- res$P

						 grad <- NULL
						 if(gradient)
								 grad <- res$gradient
						 hess <- NULL
						 if(hessian) 
								 hess <- res$hessian

						 attributes(f) <- list(C = C, A = A, P = P, gradient=grad, hessian=hess)
						 f
}
				)
}
		)
		class(result) <- "semObjective"
		result
}

#objectivelogLik <- function(gradient=FALSE, hessian=FALSE)
#{
#		result<-list(objective = function(par, object){
#								 with(object, {

#											res <- CompiledObjective(par=par, model.description=object, objective="objectivelogLik", gradient=gradient, hessian=hessian)

#		f <- res$f
#		grad <- NULL
#		if(gradient)
#				grad <- res$gradient
#		hess <- NULL
#		if(hessian)
#				hess <- res$hessian
#		attributes(f) <- list(gradient=grad, hessian=hess)
#		f
#				}
#		)
#		}
#		)
#		result
#}


logLik.objectiveFIML <- function(object, saturated=FALSE, intercept="Intercept", iterlim=1000, ...){
		logLikSaturated <- function(object, iterlim, ...){
		#objective <- function(par)
		#{

		#		res <- CompiledObjective(par=par, model.description=object, objective="objectivelogLik", gradient=FALSE, hessian=FALSE)

		#		f <- res$f
		#		f
		#}

		data <- object$data
		valid <- !is.na(data)
		valid.pattern <- apply(valid, 1, function(row) paste(row, collapse="."))
		unique.patterns <- unique(valid.pattern)
		pattern.number <- apply(outer(valid.pattern, unique.patterns, `==`), 1, which)
		valid.data.patterns <- t(sapply(strsplit(unique.patterns, "\\."), as.logical))
		n.pat <- nrow(valid.data.patterns) 
		log.2pi <- log(2*pi)
		n <- ncol(data)
		N <- nrow(data)
		C <- object$C
		tri <- lower.tri(C, diag=TRUE)

		posn.intercept <- which(rownames(C) == intercept)

		tri[posn.intercept, posn.intercept] <- FALSE
		start <- C[tri]
		opt <- options(warn=-1)
		on.exit(options(opt))
		object$posn.intercept <- posn.intercept
		object$intercept <- intercept
		object$pattern.number <- pattern.number
		object$valid.data.patterns <- valid.data.patterns
		object$tri <- tri

		#res <- nlm(objective, start, iterlim=iterlim)
		res <- CompiledSolve(model.description=object, start=start, objective="objectivelogLik",  maxiter=iterlim)

		logL <- - res$minimum/2
		C <- matrix(0, n, n)
		C[tri] <- res$estimate
		C <- C + t(C) - diag(diag(C))
		C[posn.intercept, posn.intercept] <- 1
		list(logL=logL, C=C, code=res$code)
		}
		if (saturated) {
				res <- logLikSaturated(object, iterlim=iterlim)
				if (res$code > 3) warning("nlm return code = ", res$code)
				logL <- res$logL
				attr(logL, "C") <- res$C
				return(logL)
		}
		else return(- object$criterion*object$N/2)
}

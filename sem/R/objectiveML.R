# last modified 2012-01-06 by J. Fox
# Modified for Compiled code in C/C++ by Zhenghua Nie.

objectiveML <- function(gradient=TRUE, hessian=FALSE){
	result <- list(
		objective = function(par, model.description){
			with(model.description, {

					res <- CompiledObjective(par=par, model.description=model.description, objective="objectiveML",gradient=gradient, hessian=hessian) 

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
					attributes(f) <- list(C=C, A=A, P=P, gradient=grad, hessian=hess)
					f
}
)
		}
)
if (gradient)
		result$gradient <- function(par, model.description){
				with(model.description, {

					res <- CompiledObjective(par=par, model.description=model.description, objective="objectiveML",hessian=hessian) 

				  A <- res$A
					P <- res$P
					C <- res$C
					grad <- res$gradient

					attributes(grad) <- list(C=C, A=A, P=P, gradient=grad)
					grad
				}
			)
		}
if (hessian)
		result$hessian <- function(par, model.description){
				with(model.description, {

					res <- CompiledObjective(par=par, model.description=model.description, objective="objectiveML",hessian=hessian) 

				  A <- res$A
					P <- res$P
					C <- res$C
					hess <- res$hessian

					attributes(grad) <- list(C=C, A=A, P=P, hessian=hess)
					grad
				}
			)
		}
	class(result) <- "semObjective"
	result
}

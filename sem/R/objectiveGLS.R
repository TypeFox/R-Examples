# last modified 2012-01-06 by J. Fox
# Modified for Compiled code in C/C++ by Zhenghua Nie.

objectiveGLS <- function (gradient = FALSE) 
{
    result <- list(objective = function(par, model.description) {
        with(model.description, {
						 
						 res <- CompiledObjective(par=par, model.description=model.description, objective="objectiveGLS", gradient=gradient) 
				
						 f <- res$f
				     C <- res$C
				     A <- res$A
				     P <- res$P

            attributes(f) <- list(C = C, A = A, P = P)
            f
        }
				)
    }
		)
    class(result) <- "semObjective"
    result
}

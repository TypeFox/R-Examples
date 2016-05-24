# last modified 2012-01-11 by J. Fox

objectiveML2 <- function(gradient=TRUE){
	result <- list(
		objective = function(par, model.description){
			with(model.description, {
					A <- P <- matrix(0, m, m)
					val <- ifelse (fixed, ram[,5], par[sel.free])
					A[arrows.1] <- val[one.head]
					P[arrows.2t] <- P[arrows.2] <- val[!one.head]
					I.Ainv <- solve(diag(m) - A)
					C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
					Cinv <- solve(C)
					f <- sum(diag(S %*% Cinv)) + log(det(C)) - n - logdetS
					grad <- NULL
					if (gradient){
						grad.P <- correct * t(I.Ainv) %*% t(J) %*% Cinv %*% (C - S) %*% Cinv %*% J %*% I.Ainv
						grad.A <- grad.P %*% P %*% t(I.Ainv)        
						grad <- rep(0, t)
						grad[sort(unique.free.1)] <- tapply(grad.A[arrows.1.free],ram[ram[,1]==1 & ram[,4]!=0, 4], sum)
						grad[sort(unique.free.2)] <- tapply(grad.P[arrows.2.free],ram[ram[,1]==2 & ram[,4]!=0, 4], sum)
					}
					attributes(f) <- list(C=C, A=A, P=P, gradient=grad)
					f
				}
			)
		}
	)
	if (gradient)
		result$gradient <- function(par, model.description){
			with(model.description, {
					A <- P <- matrix(0, m, m)
					val <- ifelse (fixed, ram[,5], par[sel.free])
					A[arrows.1] <- val[one.head]
					P[arrows.2t] <- P[arrows.2] <- val[!one.head]
					I.Ainv <- solve(diag(m) - A)
					C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
					Cinv <- solve(C)
					grad.P <- correct * t(I.Ainv) %*% t(J) %*% Cinv %*% (C - S) %*% Cinv %*% J %*% I.Ainv
					grad.A <- grad.P %*% P %*% t(I.Ainv)        
					grad <- rep(0, t)
					grad[sort(unique.free.1)] <- tapply(grad.A[arrows.1.free],ram[ram[,1]==1 & ram[,4]!=0, 4], sum)
					grad[sort(unique.free.2)] <- tapply(grad.P[arrows.2.free],ram[ram[,1]==2 & ram[,4]!=0, 4], sum)
					attributes(grad) <- list(C=C, A=A, P=P, gradient=grad)
					grad
				}
			)
		}
	class(result) <- "semObjective"
	result
}

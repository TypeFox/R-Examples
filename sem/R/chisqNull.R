chisqNull <- function(object){
	UseMethod("chisqNull")
}

chisqNull.objectiveML <- function(object){
	chisq <- if (!object$raw) {
			S <- object$S
			CC <- diag(diag(S))
			(object$N - 1) * 
				(sum(diag(S %*% solve(CC))) + log(det(CC)) - log(det(S)) - object$n)
		}
		else NULL
	chisq
}

chisqNull.objectiveGLS <- function(object){
	chisq <- if (!object$raw) {
			S <- object$S
			CC <- diag(diag(S))
			SS <- solve(S) %*% (S - CC)
			(object$N - 1)*0.5*sum(diag(SS %*% SS))
		}
		else NULL
	chisq
}
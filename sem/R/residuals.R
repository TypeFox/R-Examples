# last modified 2011-11-04 by J. Fox

residuals.sem <- function(object, ...){
    object$S - object$C
    }
	
standardized.residuals <- function(...){
	.Deprecated("standardizedResiduals", package="sem")
	standardizedResiduals(...)
}

standardizedResiduals <- function(object, ...){
    UseMethod("standardizedResiduals")
    }

standardizedResiduals.sem <- function(object, ...){
    res <- residuals(object)
    s <- diag(object$S)
    res/sqrt(outer(s, s))
    }

normalized.residuals <- function(...){
	.Deprecated("normalizedResiduals", package="sem")
	normalizedResiduals(...)
    }
	
normalizedResiduals <- function(object, ...){
	UseMethod("normalizedResiduals")
}
    
normalizedResiduals.objectiveML <- function(object, ...){
    res <- residuals(object)
    N <- object$N - (!object$raw)
    C <- object$C
    c <- diag(C)
    res/sqrt((outer(c,c) + C^2)/N)
    }
	
normalizedResiduals.objectiveGLS <- function(object, ...){
	res <- residuals(object)
	N <- object$N - (!object$raw)
	S <- object$S
	s <- diag(S)
	res/sqrt((outer(s,s) + S^2)/N)
}
	

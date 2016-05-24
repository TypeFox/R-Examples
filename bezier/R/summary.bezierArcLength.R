summary.bezierArcLength <- function(object, ...){
	r <- ''
	
	r <- c(r, '\nbezierArcLength Summary\n')
	r <- c(r, '\tArc length: ', format(object$arc.length), '\n')

	if(!is.null(object$slope.break)){
		if(length(object$slope.break) == 1){
			r <- c(r, '\tSlope break: ', format(object$slope.break), '; ', object$break.cause, '', '\n')
			r <- c(r, '\tNumber of iterations: ', object$n.iter, '\n')
		}else{
			r <- c(r, '\tSlope break:\n')
			for(i in 1:length(object$slope.break)) r <- c(r, '\t\t', i, ') ', format(object$slope.break[i]), '; ', object$break.cause[i], '', '\n')

			r <- c(r, '\tNumber of iterations:\n')
			for(i in 1:length(object$slope.break)) r <- c(r, '\t\t', i, ') ', object$n.iter[i], '\n')
		}
	}

	if(!is.null(object$n)){
		if(length(object$slope.break) == 1){
			r <- c(r, '\tn: ', object$n, '\n')
		}else{
			r <- c(r, '\tn:\n')
			for(i in 1:length(object$slope.break)) r <- c(r, '\t\t', i, ') ', object$n[i], '\n')
		}
	}

	class(r) <- "summary.bezierArcLength"
	r
}

print.summary.bezierArcLength <- function(x, ...) cat(x, sep='')
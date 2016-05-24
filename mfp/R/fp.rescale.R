fp.rescale <- function(x)
{
#
	if(any(x$scale[,1]>0) | is.null(x$coefficients)) 
		return(list(coefficients=x$coefficients, variance=x$var))
#
	cox <- x$family["family"] == "Cox"
# back-transformation
    x$scaled.coefficients <- x$coefficients
    if (length(x$df.final) > 1) {
		xs <- na.omit(as.vector(t(exp(x$powers*log(x$scale[,2])))))
		if(cox) A <- diag(1/xs, nrow=length(xs)) else A <- diag(c(1,1/xs))
		# log-trafo
			scales <- as.vector(t(x$scale[,2]*(!is.na(x$powers)))); scales <- scales[scales>0]
			logtr <- which(na.omit(as.vector(t(x$powers)==0)))
			if(!cox) if(length(logtr)) A[1,logtr+1] <- -log(scales[logtr])
        # equal trafos 
		logtr <- cumsum(!is.na(t(x$power)))[which(x$power[,2]==x$power[,1])*2]
			if(length(logtr)) { 
				if(cox) sel <- logtr-1 else sel <- logtr
				A[sel, sel+1] <- -log(scales[sel])/xs[sel] 
			}
			x$scaled.coefficients <- as.vector(A %*% x$coefficients); names(x$scaled.coefficients) <- names(x$coefficients)
			x$scaled.var <- A %*% x$var %*% t(A); dimnames(x$scaled.var) <- list(names(x$scaled.coefficients), names(x$scaled.coefficients))
    }
    else {
        if(x$df.final != 0) {
            p <- 1; if(x$df.final == 4) p <- 2
			xs <- exp(x$powers[1:p]*log(x$scale[2])) # power transformed scaling factors (1 for log-trafo)
			if(cox) A <- diag(1/xs, nrow=length(xs)) else A <- diag(c(1,1/xs))
			# log-trafo
			scales <- as.vector(t(x$scale[,2]*(!is.na(x$powers)))); scales <- scales[scales>0]
			logtr <- which(x$powers==0)
			if(!cox) if(length(logtr)) A[1,logtr+1] <- -log(scales[logtr])
			# equal trafos 
			logtr <- which(x$powers[,2]==x$powers[,1])*2
			
			if(length(logtr)) { 
				if(cox) sel <- logtr-1 else sel <- logtr
				A[sel, sel+1] <- -log(scales[sel])/xs[sel] 
			}
			#
			x$scaled.coefficients <- as.vector(A %*% x$coefficients); names(x$scaled.coefficients) <- names(x$coefficients)
			x$scaled.var <- A %*% x$var %*% t(A); dimnames(x$scaled.var) <- list(names(x$scaled.coefficients), names(x$scaled.coefficients))
        }
    }
    return(list(coefficients=x$scaled.coefficients, variance=x$scaled.var))
}

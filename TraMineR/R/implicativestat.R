###########################
## Compute the implicative statistic of a rule
###########################


implicativestat <- function(x, y, type="intensity", resid="standard") {
	if (!(type %in% c("intensity", "indice"))) {
		stop("type should be intensity or indice")
	}
	if (!(resid %in% c("standard", "deviance", "Freeman-Tukey", "adjusted"))) {
		stop("resid should be one of standard, deviance, Freeman-Tukey or ajusted")
	}
	x <- factor(x)
	y <- factor(y)
	xgrp <- levels(x)
	ygrp <- levels(y)
	result <- matrix(0, nrow=length(xgrp), ncol=length(ygrp))
	n <- length(x)
	rownames(result) <- xgrp
	colnames(result) <- ygrp
	for (i in 1:length(xgrp)) {
		condi <- x==xgrp[i]
		for (j in 1:length(ygrp)) {
			condj <- y==ygrp[j]
			Nnbj <- sum((condi)&!condj)
			Nexpnbj <- sum(condi)*sum(!condj)/n
			if (resid=="standard") {
				indice <- (Nnbj-Nexpnbj)/sqrt(Nexpnbj)
			}
			else if (resid=="deviance") {
				indice <- sign(Nnbj-Nexpnbj)*sqrt(abs(2*Nnbj*log(Nnbj/Nexpnbj)))
			}
			else if (resid=="Freeman-Tukey") {
				indice <- sqrt(Nnbj)+sqrt(1+Nnbj)-sqrt(4*Nexpnbj+1)
			}
			else if (resid=="adjusted") {
				indice <- (Nnbj-Nexpnbj)/sqrt(Nexpnbj*(sum(condj)/n)*(1-sum(condi)/n))
			}
						
			
			if (type=="indice") {
				result[i, j] <- indice
			}
			else {
				result[i, j] <- pnorm(-indice)
			}
		}
	}
	return(result)
}
	
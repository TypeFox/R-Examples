# modified 24 May 2009 by Adam Kramer (original by J. Fox)
# last modified 2011-08-10 by J. Fox

standardized.coefficients <- function(...){
	.Deprecated("standardizedCoefficients", package="sem")
	standardizedCoefficients(...)
}

standardizedCoefficients <- function(object, ...){
	UseMethod("standardizedCoefficients")
}

standardizedCoefficients.sem <- function (object, digits = getOption("digits"), oneheaded = TRUE, twoheaded = TRUE, ...) {
	if (!oneheaded && !twoheaded) {
		stop("No coefficients requested.")
	}
	old.digits <- options(digits = digits)
	on.exit(options(old.digits))
	P <- object$P
	A <- object$A
	t <- object$t
	par <- object$coeff
	par.posn <- object$par.posn
	IAinv <- solve(diag(nrow(A)) - A)
	C <- IAinv %*% P %*% t(IAinv)
	ram <- object$ram
	par.names <- rep(" ", nrow(ram))
	for (i in 1:t) {
		which.par <- ram[, 4] == i
		ram[which.par, 5] <- par[i]
		par.names[which.par] <- names(par)[i]
	}
	coeff <- ram[, 5]
	if (oneheaded) {
		one.head <- ram[, 1] == 1
		coeff[one.head] <- coeff[one.head] * sqrt(diag(C[ram[one.head, 
												3], ram[one.head, 3], drop = FALSE])/diag(C[ram[one.head, 
												2], ram[one.head, 2], drop = FALSE]))
	}
	if (twoheaded) {
		two.head <- ram[, 1] == 2
		coeff[two.head] <- coeff[two.head]/sqrt(diag(C[ram[two.head, 
												3], ram[two.head, 3], drop = FALSE]) * diag(C[ram[two.head, 
												2], ram[two.head, 2], drop = FALSE]))
	}
	var.names <- rownames(A)
	par.code <- paste(var.names[ram[, 2]], c("<---", "<-->")[ram[, 
							1]], var.names[ram[, 3]])
	coeff <- data.frame(par.names, coeff, par.code)
	colnames(coeff) <- c(" ", "Std. Estimate", " ")
	if (oneheaded && twoheaded) {
		coeff
	}
	else if (oneheaded) {
		coeff[one.head, ]
	}
	else {
		coeff[two.head, ]
	}
}

std.coef <- function(...){
	.Deprecated("stdCoef", package="sem")
	standardizedCoefficients(...)
}

stdCoef <- function (...){
    standardizedCoefficients(...)
    }

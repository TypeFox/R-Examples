#' @title Test actor and partner means for equality
#'
#' @description
#' Test actor and partner means for equality, using a Wald test.
#'
#' @export
#' @import lavaan
#' @param x A fSRM object-
#' @param digits Digits to which the printed results are rounded

equalMeans <- function(x, digits=3) {
	if (x$means == FALSE) stop("Please provide an fSRM object with mean structure (set `means` to TRUE)")
	eff <- parameterEstimates(x$fit, standardized=TRUE)
	vc <- vcov(x$fit)
	r <- x$roles[1]

	# actor means
	A <- eff[grepl(".means.A", eff$label) & eff$lhs != paste0(style$actor, ".", r), ]
	A.vc <- vc[A$label, A$label]
	Wald.A <- A$est %*% solve(A.vc) %*% A$est

	# partner means
	P <- eff[grepl(".means.P", eff$label) & eff$lhs != paste0(style$partner, ".", r), ]
	P.vc <- vc[P$label, P$label]
	Wald.P <- P$est %*% solve(P.vc) %*% P$est
	
	# relationship means
	if (length(x$roles) == 3) {
		R <- eff[grepl(".means.R", eff$label), ][1, ]
		R.vc <- vc[R$label, R$label]
		Wald.R <- R$est %*% solve(R.vc) %*% R$est
		df.R <- 1
	}
	if (length(x$roles) == 4) {
		R <- eff[grepl(".means.R", eff$label), ][c(7, 8, 5, 2, 6), ]
		R.vc <- vc[R$label, R$label]
		Wald.R <- R$est %*% solve(R.vc) %*% R$est
		df.R <- 5
	}
	
	res <- data.frame(
		Wald=c(Wald.A, Wald.P, Wald.R), 
		df = c(rep(length(x$roles)-1, 2), df.R), 
		p.value=c(1-pchisq(Wald.A, length(x$roles)-1), 1-pchisq(Wald.P, length(x$roles)-1), 1-pchisq(Wald.R, df.R)))
	rownames(res) <- c("H0: Equal actor means", "H0: Equal partner means", "H0: Equal relationship means")
	
	class(res) <- "eqM"
	return(res)
}

#' @export
#' @method print eqM
print.eqM <- function(x, ...) {
	x$sig <- p2star(x$p.value)
	x$p.value <- p(x$p.value)
	x$Wald <- round(x$Wald, 3)
	class(x) <- "data.frame"
	print(x)
}
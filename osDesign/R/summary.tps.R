summary.tps <-
function(object, ...)
{
	# produces summary from an object of the class "tps"
	coef <- object$coef
	method <- object$method
	if(method == "WL")
	{
		see <- sqrt(diag(object$cove))
		tee <- coef/see
		pee <- 2 * (1 - pnorm(abs(tee)))
		coefficients <- matrix(0, nrow = length(coef), ncol = 4)
		dimnames(coefficients) <- list(names(coef), c("Value", "Emp SE", "Emp t", "Emp p"))
		coefficients[, 1] <- coef
		coefficients[, 2] <- see
		coefficients[, 3] <- tee
		coefficients[, 4] <- pee
	}
	else
	{
		se  <- sqrt(diag(object$covm))
		see <- sqrt(diag(object$cove))
		te  <- coef/se
		tee <- coef/see
		pe  <- 2 * (1 - pnorm(abs(te)))
		pee <- 2 * (1 - pnorm(abs(tee)))
		coefficients <- matrix(0, nrow = length(coef), ncol = 7)
		dimnames(coefficients) <- list(names(coef), c("Value", "Mod SE", "Mod t", "Mod p", "Emp SE", "Emp t", "Emp p"))
		coefficients[, 1] <- coef
		coefficients[, 2] <- se
		coefficients[, 3] <- te
		coefficients[, 4] <- pe
		coefficients[, 5] <- see
		coefficients[, 6] <- tee
		coefficients[, 7] <- pee
	}
	structure(list(coefficients = coefficients), class = "summary.tps")
}

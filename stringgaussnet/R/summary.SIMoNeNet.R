summary.SIMoNeNet <-
function (object, ...)
{
	cat("Thetas summary:\n")
	print(summary(object$Edges$Theta))
	cat("\nRhos summary:\n")
	print(summary(object$Edges$Rho))
	cat("\nAbsolute rhos summary:\n")
	print(summary(abs(object$Edges$Rho)))
	cat("\nP-values summary:\n")
	print(summary(object$Edges$P.Value))
}

print.bigRR <-
function(x, ...) {
	Obj <- x
	opt <- options()
	dgt <- opt$digits
	cat("Call: \n")
	print(Obj$Call)
	cat("\n")
	cat("Fixed effects estimates:\n")
	print(Obj$beta)
	cat("\n")
	cat("Shrinkage estimates of the random effects (Genetic effects of", length(as.numeric(Obj$u)), "markers):\n")
	print(summary(as.numeric(Obj$u)))
	cat("\n")
	cat("Random effects (Genetic) variance component estimate:\n")
	print(Obj$lambda)
	cat("\n")
	cat("Residual variance component estimate:\n")
	print(Obj$phi)
	cat("\n")
}


summary.vscc <- function(object, ...){
	x <- object
    cat("---------- Summary for VSCC ----------", "\n\n")
	cat("         ----   RESULTS   ----    ", "\n")
	cat("# Vars:     ", ncol(x$top), "\n")
	cat("Relation:   ", x$chosen, "\n")
	cat("BIC:        ", x$bestmod$bic, "\n")
	cat("Model:      ", x$bestmod$model, "\n")
	cat("Family:     ", x$family, "\n")
	cat("# Groups:   ", x$bestmod$G, "\n")
}

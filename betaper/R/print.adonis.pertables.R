`print.adonis.pertables` <-
function(x,...) {
m <- x$raw
n <- x$simulation$R2.quant
t <- x$simulation$p.quant
cat("Permutation tests for multivariate analysis of variance using distance matrices", "\n")
print(m)
cat("\n\n")
cat("Confidence intervals of R-squared under different taxonomic scenarios", "\n\n")
print(n)
cat("\n\n")
cat("Confidence intervals of p-values under different taxonomic scenarios", "\n\n")
print(t)
}


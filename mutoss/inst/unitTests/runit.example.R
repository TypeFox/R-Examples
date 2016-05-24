test.bonferroni <- function() {
	result <- bonferroni(c(0.1,0.2,0.3,0.4))
	# Check for correct results:
	checkEquals(result$adjPValues, c(0.4, 0.8, 1.0, 1.0))
	# We can also check for exceptions:
	# checkException(bonferroni(c(-0.1,0.2,1.3)))
	# But unfortunately there is no check for valid pvalues yet.
}
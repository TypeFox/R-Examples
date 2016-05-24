print.summary.BiDimRegression <-
function(x, ...) 
{
	# -- short version of the statistics
	# Overall analysis
	overallStats=matrix(c(
            x$euclidean.r, x$euclidean.rsqr, x$euclidean.fValue, x$euclidean.df1, x$euclidean.df2, x$euclidean.pValue, 
		x$affine.r, x$affine.rsqr, x$affine.fValue, x$affine.df1, x$affine.df2, x$affine.pValue),
		nrow=2, byrow = TRUE)
      colnames(overallStats) <- c("r", "r-sqr", "F-value", "df1", "df2", "p-value")
	rownames(overallStats) <- c("Euclidean", "Affine")

	cat(sprintf('\n--- summary statistics from bidimensional regressions ---\n'))
	printCoefmat(overallStats, P.values=TRUE, has.Pvalue=TRUE, digits=3, tst.ind=4, signif.stars=TRUE)

}

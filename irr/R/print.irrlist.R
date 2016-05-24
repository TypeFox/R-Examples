"print.irrlist" <-
function(x, ...)
{
  cat(" ", x$method, "\n\n",sep="")
	cat(paste(" Subjects =", x$subjects, "\n"))
	cat(paste("   Raters =", x$raters, "\n"))
	results <- paste(formatC(x$irr.name, width=9, flag="+"), "=", format(x$value, digits=3),"\n")
	cat(results)
	if (!is.null(x$statistic)) {
  	statistic <- paste(formatC(x$stat.name, width=9, flag="+"), "=", format(x$statistic, digits=3), "\n")
  	cat("\n", statistic, sep="")
  	cat(paste("  p-value =", format(x$p.value, digits=3), "\n"))
	}
  if (!is.null(x$detail)) {
    cat("\n")
    if (is.table(x$detail)) print(x$detail)
    else cat(x$detail,"\n")
  }
  if (!is.null(x$error)) cat("\n ", x$error, "\n", sep="")
}


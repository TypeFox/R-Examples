"print.icclist" <-
function(x, ...)
{
  icc.title <- ifelse(x$unit=="single", "Single Score Intraclass Correlation", "Average Score Intraclass Correlation")
	cat(paste(" ",icc.title,"\n\n",sep=""))
	cat(paste("   Model:", x$model, "\n"))
	cat(paste("   Type :", x$type, "\n\n"))
	cat(paste("   Subjects =", x$subjects, "\n"))
	cat(paste("     Raters =", x$raters, "\n"))
	results <- paste(formatC(x$icc.name, width=11, flag="+"), "=", format(x$value, digits=3))
	cat(results)
  cat("\n\n F-Test, H0: r0 =",x$r0, "; H1: r0 >",x$r0,"\n")
	Ftest <- paste(formatC(paste("F(",x$df1,",",format(x$df2, digits=3),")",sep=""), width=11, flag="+"), "=",
                                              format(x$Fvalue, digits=3),
	               ", p =", format(x$p.value, digits=3), "\n\n")
  cat(Ftest)
	cat(" ", round(x$conf.level*100,digits=1), "%-Confidence Interval for ICC Population Values:\n", sep="")
	cat(paste("  ", round(x$lbound, digits=3), " < ICC < ", round(x$ubound, digits=3), "\n", sep=""))
}


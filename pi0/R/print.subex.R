`print.subex` <-
function(x,...)
{
print(x$extrp.fit)
cat("\np-values summary",fill=TRUE)
print(quantile(x$pvalues,c(0:5/100,2:10/20,1)),digits=3)
cat("\nq-values summary",fill=TRUE)
print(quantile(x$qvalues,c(0:5/100,2:10/20,1)),digits=3)
invisible(NULL)
}


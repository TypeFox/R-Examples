print.perm.oneway.anova <-
function(x, ...)
{
hist(x$Perm.values,breaks=20,xlab="permutation F-value",
 main="Histogram of permutation F-values")
abline(v=x$Observed,col="2")
leg.text <- expression(Observed)
legend("topright",leg.text,col=2,lwd=2,cex=.6)
cat("\n\n",x$Header,"\n\n")
cat("SUMMARY STATISTICS\n\n")
u <- data.frame(x$Levels,n=x$n,Mean=x$Mean,SD=x$SD)
names(u)[1]=x$Factor
print(u,row.names=FALSE)
cat("\n")
cat("HYPOTHESIS TEST\n\n")
print(data.frame(Response=x$Response,Factor=x$Factor,Trim=x$Trim,Statistic=x$Statistic,
 Observed=x$Observed,P.value=x$P.value),row.names=FALSE)
cat("\n\n")
}

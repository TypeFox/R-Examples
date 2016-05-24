print.perm.two.var <-
function(x, ...)
{
hist(x$Perm.values,breaks=20,xlab=paste("permutation",x$Statistic),
  main=paste("Histogram of permutation ",x$Statistic,"s",sep=""))
abline(v=x$Observed,col="2")
leg.text <- expression(Observed)
legend("topright",leg.text,col=2,lwd=2,cex=.6)
cat("\n\n",x$Header,"\n\n")
print(data.frame(SUMMARY="STATISTICS",Variable.1=x$Variable.1,
 Variable.2=x$Variable.2,n=x$n,Statistic=x$Statistic,
  Observed=x$Observed),row.names=F)
cat("\n")
print(data.frame(HYPOTHESIS="TEST",Null=x$Null,Alternative=x$Alternative,
  P.value=x$P.value),row.names=F)
cat("\n\n")
}

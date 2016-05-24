print.perm.paired.loc <-
function(x, ...)
{
hist(x$Perm.values,breaks=20,xlab=paste("permutation",x$Statistic),
  main=paste("Histogram of permutation ",x$Statistic,"s",sep=""))
abline(v=x$Observed,col="2")
leg.text <- expression(Observed)
legend("topright",leg.text,col=2,lwd=2,cex=.6)
cat("\n\n",x$Header,"\n\n")
if (!is.null(x$Variable))
 print(data.frame(SUMMARY="STATISTICS",Variable=x$Variable,
  Pop.1=x$Pop.1,Pop.2=x$Pop.2,n=x$n,Statistic=x$Statistic,
  Observed=x$Observed),row.names=FALSE)
else
 print(data.frame(SUMMARY="STATISTICS",Pop.1=x$Pop.1,Pop.2=x$Pop.2,n=x$n,
 Statistic=x$Statistic,Observed=x$Observed),row.names=FALSE)
cat("\n")
print(data.frame(HYPOTHESIS="TEST",Null=x$Null,Alternative=x$Alternative,
 P.value=x$P.value),row.names=FALSE)
cat("\n\n")
}

print.rmpobject <- function(x, ...) 
{
cat(x$mixing_type,x$name,"\n")
cat(x$mcmc$nrep, "iterations completed in ", as.double(x$mcmc$time), attr(x$mcmc$time, "units"), "\n")
cat("Posterior mean number of cluster ", as.double(x$clustering$groups), "\n")
}
#----------------------------------------
summary.rmpobject <- function(object, ...) 
{
cat("------------------------------\n",
	object$mixing_type, object$name, "\n------------------------------\n",
   "Posterior mean number of cluster ", object$clustering$groups,
	"\nPosterior pmf and 95% pointwise credible interval: \n", sep="")
print(data.frame(lower = round(object$pmf$lower.95,4), pmf=round(object$pmf$post.pmf,4), upper = round(object$pmf$upper.95,4), row.names=object$pmf$domain))
cat("-----------------------\n")
}
#---
plot.rmpobject <- function(x, y=NULL, main="", col=2, nb=x$mcmc$nb, ...) 
{
post.pmf<-apply(x$mcmc.chains$pmf[(nb+1):x$mcmc$nrep,],2,mean)[x$pmf$lb:x$pmf$ub+1]
sorted.p<-apply(x$mcmc.chains$pmf[(nb+1):x$mcmc$nrep,],2,sort)
pmf2.5 <-sorted.p[(0.025*(x$mcmc$nrep-nb)),][x$pmf$lb:x$pmf$ub+1]
pmf97.5 <-sorted.p[(0.975*(x$mcmc$nrep-nb)),][x$pmf$lb:x$pmf$ub+1]
plot(x$pmf$empirical~x$pmf$domain, lty=2, ty='h', main=main, ylab="pmf", xlab="estimated pmf and empirical pmf (dashed)")
points(x$pmf$post.pmf~c(0.2+x$pmf$domain), col=col, ty='h')
}


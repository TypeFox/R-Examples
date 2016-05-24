plot.path = function(x, ...)
{
	par(mfrow=c(2,1), mar=c(4,5,2,2), cex.axis=0.75)
	colors = 3:(nrow(x$beta)+2)
	ub = max(x$breaks)*1.08

	plot(1, type="n", xlab=expression(lambda), ylab=expression(beta), xaxt="n", 
	     ylim=range(x$beta), xlim=c(0, ub))
    axis(1, at=round(x$breaks,2), labels=round(x$breaks, 2))
	for (i in 1:nrow(x$beta)) 
		points(c(ub,x$breaks), c(x$beta[i,1],x$beta[i,]), type="l", col=colors[i], lwd=2)
    abline(0,0)
		
	plot(1, type="n", xlab=expression(lambda), ylab="Score", xaxt="n",
	     ylim=c(-max(x$score), max(x$score)), xlim=c(0, ub))
    axis(1, at=round(x$breaks,2), labels=round(x$breaks, 2))
    points(c(0, ub), c(0, ub), col=2, type="l", lwd=3)
    points(c(0, ub), c(0, -ub), col=2, type="l", lwd=3)
    abline(0,0)

	for (i in 1:nrow(x$score))
		points(c(ub,x$breaks), c(x$score[i,1],x$score[i,]), type="l", col=colors[i], lwd=2)
}

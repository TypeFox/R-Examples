plot.SOGLresult<-function(x, which, ...) {
if (!(which %in% c("alpha selection", "density", "empirical CI")) ) stop("which must be one of: alpha selection, density, empirical CI")
if (which == "alpha selection")
{
 plot(x$alpha.considered,x$pAUC, type="l", ylab="pAUC score", xlab=expression(alpha),log="x")
 abline(v=x$alpha.selected, col="red")
 mtext(expression(alpha[sel]), side = 3, at = x$alpha.selected)
 }
if (which == "density")
{
  ran<-range(x$random, x$subsample, x$score)
  drand<-density(x$random)
  dsub<-density(x$subsample)
  ran2<-range(drand$y, dsub$y)*1.2
  
  plot(drand, main="", xlab="Similarity score", ylab="Density", xlim=ran , ylim=ran2)
  lines(dsub, col = 2)
  abline(v=x$score, col = 2)
  mtext("observed", side= 3, at= x$score)
  text(x=max(ran), y=ran2[2], adj=c(1,1), labels=paste(paste("p-value: ", round(x$significance,3)), 
                            paste("alpha: ", round(x$alpha.selected,3)),sep="\n") )
  rug(x$subsample, col=2)
  rug(x$random)                     
  legend(x="topleft", legend=c("subsampled","random"), col = c(2,1), lty = 1, bty="n")
}      
if (which == "empirical CI")
{
 N<-length(x$common.gene$top) 
 matplot(rbind(x$emp.ci$top[1:N,], flip(x$emp.ci$bottom[1:N,])), type="n", xaxt="n",ylab="Overlaping genes")
 for (i in 1:N) lines(x=c(i,i), y=x$emp.ci$top[i,c(1,3)], col="orange")
 lines( x$emp.ci$top[1:N,2], col="orange", type="S")
 lines(x$common.gene$top[1:N], col="black", type="S")  
bottom<-flip(x$emp.ci$bottom[1:N,])
 for (i in ((1:N)+N)) lines(x=c(i,i), y=bottom[i-N,c(1,3)], col="orange")
 lines( x=c(1:N+N), y=bottom[1:N,2], col="orange", type="S")
 lines( x=c(1:N+N), y=rev(x$common.gene$bottom[1:N]), col="black", type="S")

abline(v=N, lty=2)
legend(x= "topleft", legend=c("Observed", "Expected"), col=c("black","orange"), lty=1, bty="n")  
axStep <- floor(N/45) * 10
axis(side = 1, at = c(axStep * (0:4), 2 * N - axStep * (0:4)), labels = c(1, axStep * (1:4), 1, axStep * (1:4)))
mtext("Top ranks", side = 1, at = 0.5 * N, line = 2.5)
mtext("Bottom ranks", side = 1, at = 1.5 * N, line = 2.5) 
}
}




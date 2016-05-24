roc.curves.plot <-
function(sim1.curve,sim2.curve,mod1,mod2) {

# draw the original curve
plot(sim1.curve,col="blue",type="o",pch=16,lwd=2,yaxs="i",xaxs="i")
plot(sim2.curve,col="red",type="o",pch=20, lwd=2,lty=3, add=TRUE)
abline(0,1,col="gray",lty=2)
legend(0.7,0.7,c(mod1,mod2),pch=c(16,20),col=c("blue","red"),lwd=2,lty=c(1,3), bty = "n",box.lty=0)
mtext(paste("Empirical ROC curves"),line=3)
}

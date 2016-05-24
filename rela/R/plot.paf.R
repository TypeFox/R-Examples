plot.paf <- function(x, ...) {

# Plotting reproduced correlation differences

dev.new()
if (is.integer(nrow(x$Residuals)/2)==FALSE) {par.lv <- (nrow(x$Residuals)/2)+.5 }
else if (is.integer(nrow(x$Residuals)/2)==TRUE) {par.lv <- nrow(x$Residuals)/2 }
par(mfcol=c(par.lv,2), mar=c(2, 2, 3, 1), oma=c(2, 1, 4, 1))
for (i in 1:ncol(x$Residuals) ) {
yl <- round(min(x$Residuals[,i]),digits=2)
yh <- round(max(x$Residuals[,i]),digits=2)
plot(round(x$Residuals[,i], digits=2), type="b", pch=19, lty=1, lwd=1, col="dark green",axes = FALSE, frame.plot = TRUE,
ylim=c(yl, yh), ylab="Difference", xlab=NA, main=colnames(x$Residuals)[i] )
axis(1, at=c(1:nrow(x$Residuals)), labels=c(rownames(x$Residuals)) )
axis(2, at=c(round(min(x$Residuals[,i]),digits=2),0,round(max(x$Residuals[,i]),digits=2)) ) 
abline(h=0, lty=2, lwd=.5, col="black") }
mtext("Correlation Residuals", outer=TRUE, line=1, cex=1.2)


# Plotting communalities

dev.new()
yclma <- round(max(x$Communalities),digits=2)
yclmi <- round(min(x$Communalities),digits=2)
plot(round(x$Communalities[,1], digits=2), type="b", pch=19, lty=1, lwd=2, col="dark green", 
axes = FALSE, frame.plot = TRUE, ylim=c(yclmi, yclma), ylab="Communalities", xlab="Items", 
main="Item Communalities" )
points(round(x$Communalities[,2], digits=2), type="b", pch=19, lty=2, lwd=2, col="dark orange")
axis(1, at=c(1:nrow(x$Communalities)), labels=c(rownames(x$Communalities)) )
axis(2, at=round(c(x$Communalities[,1],x$Communalities[,2]),digits=2) )
legend(locator(1), c("Initial", "Final"), lty=c(1,2), lwd=c(2,2), col = c("dark green", "dark orange"), box.lty=0)


# Plotting MSAs

dev.new()
ylma <- round(max(x$MSA),digits=2)
ylmi <- round(min(x$MSA),digits=2)
plot(round(x$MSA, digits=2), type="b", pch=19, lty=1, lwd=2, col="dark blue",axes = FALSE, frame.plot = TRUE,
ylim=c(ylmi, ylma), ylab="MSA", xlab="Item", main="Measures of Sampling Adequacy" )
axis(1, at=c(1:nrow(x$MSA)), labels=c(rownames(x$MSA)) )
axis(2, at=round(x$MSA,digits=2) )


# Plotting eigenvalue scree plot

dev.new()
plot(x$Eigenvalues[,1], pch=19, type="o", ylab="Eigenvalues", xlab="Components", main="Principal Axis Factoring Scree Plot")
abline(h=as.numeric(as.character(x$call[3])), lty=2, lwd=2, col="red")


}
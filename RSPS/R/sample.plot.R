sample.plot <-
function(fit,ylim=c(-0.05,1.05),cutoff=0.8)
{
Mean.Null <- Type.I.Error <- Effect.Size <- NULL
fit <- fit
if (dim(fit)[2]==7) {dat.temp <- subset(fit,Mean.Null==min(Mean.Null),Type.I.Error=min(Type.I.Error));main="Power Curve (Poisson Dist.)"}
if (dim(fit)[2]==8) {dat.temp <- subset(fit,Mean.Null==min(Mean.Null),Type.I.Error=min(Type.I.Error));main="Power Curve (Neg. Bin. Dist.)"}

key <- list(columns = 3,points=FALSE,lines=TRUE,lty = 1:4)
plot <- xyplot(Power.est ~ N.est,groups = Effect.Size,data= dat.temp,lwd=2,type="l",auto.key=key,xlab="Sample Size",ylab="Power",main=main, lty = 1:length(unique(dat.temp$Effect.Size)),panel = function(...) { panel.abline(h = cutoff,col="gray",lwd=2,lty = 2);panel.xyplot(...)})
plot
}

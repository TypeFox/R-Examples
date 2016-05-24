plotprm <-
function(prmobj,y,...)
{
# Generate plots for predicted values and residuals for prm (robust PLS)
#

par(mfrow=c(1,2))
optcomp=prmobj$optcomp
plot(y,prmobj$pred[,optcomp],xlab="Measured y",ylab="Predicted y",
        cex.lab=1.2,cex=0.7,pch=3,col=1,...)
abline(c(0,1),lty=1)
resid=y-prmobj$pred[,optcomp]
ylimits=max(abs(resid))
ylimits=sort(c(-ylimits,ylimits))
plot(prmobj$pred[,optcomp],resid,
        xlab="Predicted y",ylab="Residuals",cex.lab=1.2,
        cex=0.7,pch=3,col=1,ylim=ylimits)
abline(h=0)

invisible()
}


"graph.mod" <-
function(ssmod, x, y, data, title="",xlab="Centered X",ylab="Raw Y", ylimit=1.5, ...)
{
attach(data)
mcx <- meanCenter(x)
yhi <- mean(y, na.rm=TRUE) + ylimit*sd(y, na.rm=TRUE)
ylo <- mean(y, na.rm=TRUE) - ylimit*sd(y, na.rm=TRUE)
plot(mcx,y,main=title, xlab=xlab, ylab=ylab, ylim=c(ylo, yhi))
abline(ssmod[1,1],ssmod[1,2],col=4)
abline(ssmod[3,1],ssmod[3,2],col=2)
abline(ssmod[2,1],ssmod[2,2],lty=2, col=1)
legend (locator(1), c("zHigh", "zLow"), lty=1, col=c(4,2))
detach(data)
}


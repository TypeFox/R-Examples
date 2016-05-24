plot.trace2 <-
function(x,i1,i2,j,xlab="",ylab=paste("outcome age effect at j=",j),...)
{plot(seq(i1,i2),x[i1:i2,j],type="l",xlab=xlab,ylab=ylab)}

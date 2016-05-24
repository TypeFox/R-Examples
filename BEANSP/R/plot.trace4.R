plot.trace4 <-
function(x,i1,i2,j,xlab="",ylab=paste("variance effect at j=",j),...)
{plot(seq(i1,i2),x[i1:i2,j],type="l",xlab=xlab,ylab=ylab)}

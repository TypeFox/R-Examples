plot.q <-
function(x,n,jj,xlab="nest age",ylab="individual age-specific nest failure rate",...)
{plot(seq(1,jj),x[n,1:jj],type="l",xlab=xlab,ylab=ylab)}

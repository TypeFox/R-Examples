plot.sr <-
function(x,n,jj,xlab="nest age",ylab="individual age-specific nest survival rate",...)
{plot(seq(1,jj),x[n,1:jj],type="l",xlab=xlab,ylab=ylab)}

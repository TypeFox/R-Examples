emplambdaD.fun <-
function(posE, t,lint=NULL, nint=NULL,plotEmp=TRUE, inddat=NULL, 
tit='',scax=NULL, scay=NULL)
{

if (!is.null(nint)& !is.null(lint))
stop('Error: only one of nint and lint must be provided')
if (is.null(nint)& is.null(lint))
stop('Error: one of nint and lint must be provided')
n<-length(t)
if (is.null(lint)) lint<-ceiling(n/nint)
if (is.null(nint)) nint<-ceiling(n/(lint))
if (is.null(inddat)) inddat<-rep(1,n)


indice<-rep(0,n)
indice[posE]<-1
indice[inddat==0]<-NA
int<-floor(c(0:(n-1))/lint)
emplambda<-tapply(indice, INDEX=int, FUN=mean, na.rm = TRUE)
ultlint<-n-(nint-1)*lint
pmedio1<-floor(lint/2)
pmedio<-pmedio1+c(0,cumsum(rep(lint,(nint-2))) )
pmedio<-c(pmedio, pmedio[length(pmedio)]+ceiling(ultlint/2))

cat(fill=T)
cat('Number of intervals: ', nint, ' of length: ',lint, fill=T)
if(lint*nint!=n) cat(' except the las one of length ',ultlint, fill=T)
cat(fill=T)
if (plotEmp!=FALSE)
{
if (is.null(scay)==TRUE) scay<-c(min(emplambda, na.rm=TRUE), max(emplambda, na.rm=TRUE))
if (is.null(scax)==TRUE) scax<-c(min(t, na.rm=TRUE), max(t, na.rm=TRUE))
plot(t[pmedio],emplambda,cex=0.5,xlab='time',ylab='empirical occurrence rate',
 type='s', ylim=scay, xlim=scax)
title(sub=paste('Rates calculated in ',nint,' disjoint intervals of length ',lint,sep=''), cex=0.5)
mtext(paste(tit, sep=''), outer = TRUE, line = -2,cex=1)
}

return(list(emplambda=emplambda,nint=nint, lint=lint, posE=posE, t=t))
}

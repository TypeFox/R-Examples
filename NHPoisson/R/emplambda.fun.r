emplambda.fun <-
function(posE, t, lint, plotEmp=TRUE,  inddat=NULL, tit='', scax=NULL, scay=NULL)
{
n<-length(t)
indice<-rep(0,n)
indice[posE]<-1
if (is.null(inddat)) inddat<-rep(1,n)
iini<-1
ifin<-lint

posmed<-floor(lint/2)+1
emplambda<-NULL
emplambda[1:(posmed-1)]<-NA
emplambda[posmed]<-sum(indice[iini:ifin])/sum(inddat[iini:ifin], na.rm=TRUE)
j<-posmed+1
while((ifin < n))
{
iini<-iini+1
ifin<-ifin+1
emplambda[j]<-sum(indice[iini:ifin])/sum(inddat[iini:ifin], na.rm=TRUE)
j<-j+1
}
emplambda[j:n]<-NA

if (plotEmp!=FALSE)
{
if (is.null(scay)==TRUE) scay<-c(min(emplambda, na.rm=TRUE), max(emplambda, na.rm=TRUE))
if (is.null(scax)==TRUE) scax<-c(min(t, na.rm=TRUE), max(t, na.rm=TRUE))
plot(t,emplambda,cex=0.5,xlab='time',ylab='empirical occurrence rate',
 type='s', ylim=scay, xlim=scax)
title(sub=paste('Rates  calculated in periods of length ',lint,' updated in each time unit',sep=''), cex=0.5)
mtext(paste(tit, sep=''), outer = TRUE, line = -2,cex=1)
}

return(list(emplambda=emplambda,lint=lint,posE=posE, t=t))
}

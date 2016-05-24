graphrate.fun <-
function(objres=NULL, fittedlambda=NULL, emplambda=NULL, t=NULL,
lint=NULL,typeI='Disjoint', tit='',scax=NULL,scay=NULL, xlegend='topleft',histWgraph=TRUE)
{
if (is.null(objres)&(is.null(fittedlambda)|is.null(t)|is.null(emplambda)|is.null(lint)))
stop ('Argument objres or vector of arguments (fittedlambda, emplambda, t) must be specified')

if (is.null(objres)=='FALSE')
{
fittedlambda<-objres$fittedlambda
emplambda<-objres$emplambda
typeI<-objres$typeI
      t<-objres$mlePP@t
if (typeI=='Disjoint') t<-t[objres$pm]
lint<-objres$lint
}

if (histWgraph==TRUE)	dev.new(record=TRUE)

par(mfrow=c(1,1))

if (is.null(scay)==TRUE)
{
yminn<-min(emplambda,fittedlambda,na.rm=TRUE)
ymaxx<-max(emplambda,fittedlambda,na.rm=TRUE)
scay<-c(yminn,ymaxx)
}
if (is.null(scax)==TRUE) scax<-c(min(t, na.rm=TRUE), max(t, na.rm=TRUE))


plot(t, emplambda, cex = 0.5, xlab = "time", ylab = "empirical and fitted  occurrence rates",  
ylim=scay,xlim=scax, type = "l",lty=1)
lines(t, fittedlambda, col='red')

legend(x=xlegend, legend=c('Empirical rate', 'Fitted rate'), col=c('black', 'red'), lty=c(1,1), cex=0.8)
title(sub=paste('Rates calculated in', typeI, 'intervals of length', lint, sep=' '), cex=0.7)
mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)
return(NULL)

}

histprod<-function(donnee,firstvar,lastvar=ncol(donnee),numr = 2,numc = 2,adjust=1) {

nbdesc <- lastvar-firstvar+1

xquant<-donnee[,firstvar:lastvar]
xrange<-max(apply(xquant,2,max,na.rm=TRUE))
#yrange<-max(hist(donnee[,firstvar],plot=FALSE)$density)
yrange<-max(density(donnee[,firstvar], na.rm = TRUE,adjust=adjust)$y)

for (i in 2:nbdesc){
#yrangeinter<-max(hist(donnee[,i+firstvar-1],plot=FALSE)$density)
yrangeinter<-max(density(donnee[,i+firstvar-1], na.rm = TRUE,adjust=adjust)$y)
yrange<-max(yrange,yrangeinter)
           }


mult <- nbdesc %/% (numr*numc)
if (nbdesc==(nbdesc %/% (numr*numc))*(numr*numc)) mult=mult-1
for (m in 0:mult) {
    par(mfrow = c(numr,numc))
    for (nbd in 1:(numr*numc)) {
          nb <- (m*(numr*numc)+nbd)
          if (nb <= nbdesc)       {

hist(donnee[,nb+firstvar-1],col=grey(0.9),border = grey(0.8),xlab=names(donnee[nb+firstvar-1]),main = paste("Histogram of" , names(donnee[nb+firstvar-1])),xlim=c(0,xrange),ylim=c(0,yrange),proba=TRUE)
step <- seq(from = 0, to = xrange, length = 100)
lines(step, dnorm(step, mean(donnee[,nb+firstvar-1], na.rm = TRUE),sd(donnee[,nb+firstvar-1], na.rm = TRUE)),lty=2)
lines(density(donnee[,nb+firstvar-1], na.rm = TRUE,adjust=adjust), lwd = 1,col="red")
                              } 
                   }
if (m < mult) dev.new()
          }  #for (m in 0:mult) {
                                  }

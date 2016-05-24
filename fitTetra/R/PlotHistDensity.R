PlotHistDensity <-
function(intens, mixres, xlim=c(0,1), trafo="asr",
                         maintitle=NULL, subtitle=NULL, xlabel=NULL, xaxis="s",
                         nbreaks=40,frequ=T,fitcol="darkgreen") {
  nbars <- ceiling(nbreaks/(max(intens,na.rm=T)-min(intens,na.rm=T)))
  h <- hist(intens, breaks=seq(0,1,by=1/nbars), plot=FALSE )
  if (frequ==TRUE) { maxy <- max(h$counts) }
    else { maxy <- max(h$density) }
  h$ylim <- c(0,1.1*maxy) #so ylim will be part of the return value
  barplot (h$counts,
        width=(h$breaks[length(h$breaks)]-h$breaks[1])/(length(h$breaks)-1),
        space=0,xlim=xlim,ylim=h$ylim,col="white",
        main=maintitle,sub=subtitle,xlab=xlabel,xaxt=xaxis,ylab="frequency")
  if (xaxis=="s") axis(1,labels=TRUE, at=seq(0,1, by=0.2))
  area <- 1
  if (frequ == TRUE) {
    classwidth<-h$breaks[2]-h$breaks[1]
    area <- length(intens[!is.na(intens)])*classwidth
  }

  mu <- mixres$psi$mu; sigma <- mixres$psi$sigma; prior <- mixres$psi$p
  ng <- length(mu)

  r <- seq(from=xlim[1], to=xlim[2], by=((xlim[2]-xlim[1])/300))
  if (trafo==1)
    for (i in 1:ng) lines(cbind(r,area*prior[i]*dnorm(r,mean=mu[i], sd=sigma[i])), lwd=1, col=fitcol )

  if (trafo=="asr") {
     for (i in 1:ng) {
      area.i= integrate(dnormasr, lower=0.001, upper=0.999 , mu=mu[i], sigma=sigma[i], subdivisions=300)$value  # the area is not 1 anymore, correct for this
      lines(cbind(r,area*prior[i]*dnorm(asin(sqrt(r)), mean=mu[i], sd=sigma[i])/area.i), lwd=1, col=fitcol )
    }
    muback<-sin(mu)
    muback<-muback*muback
    points(rep(0,ng)~muback, pch=17, col=fitcol)
  }
  else if (trafo=="logit") {
     for (i in 1:ng) {
      area.i= integrate(dnormlogit, lower=0.001, upper=0.999 , mu=mu[i], sigma=sigma[i], subdivisions=300)$value  # the area is not 1 anymore, correct for this
      lines(cbind(r,area*prior[i]*dnorm(log(r/(1-r)), mean=mu[i], sd=sigma[i])/area.i), lwd=1, col=fitcol )
    }
  }
  h # return the histogram data
}

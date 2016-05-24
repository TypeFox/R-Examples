## Function for plotting the SIR curves.  Uses list object returned by SimulateEpidemic()
PlotEpi <- function(simoutput,add=FALSE,showd=TRUE,showv=TRUE, ...)
  {
    mycols <- c(rgb(55/255,126/255,184/255),rgb(228/255,26/255,28/255),
                rgb(77/255,175/255,74/255),rgb(152/255,78/255,163/255),rgb(255/255,127/255,0/255))
    if(add)
      {
        mylty <- 2
        lines(simoutput$S,lty=mylty,col=mycols[1],lwd=2)
      }
    else
      {
        mylty <- 1
        plot(simoutput$S,type='l',ylab="individuals",xlab="time",
             col=mycols[1],lty=mylty,lwd=2,bty='n',
             ylim=c(0,max(simoutput$S)), ...)
      }
    lines(simoutput$I,lty=mylty,col=mycols[2],lwd=2)
    lines(simoutput$R,lty=mylty,col=mycols[3],lwd=2)
    if(showd)
      {
        lines(simoutput$D,lty=mylty,col=mycols[4],lwd=2)
      }
    if(showv)
      {
        lines(simoutput$V,lty=mylty,col=mycols[5],lwd=2)
      }
    if(showd & showv)
      {
        legend("topright",c("S","I","R","D","V"),fill=mycols)
      }
    if(showd & !showv)
      {
        legend("topright",c("S","I","R","D"),fill=mycols[1:4])
      }
    if(showv & !showd)
      {
        legend("topright",c("S","I","R","V"),fill=mycols[c(1:3,5)])
      }
    if(!showd & !showv)
      {
        legend("topright",c("S","I","R"),fill=mycols[1:3])
      }
  }

## function for plotting the cost trajectories associated with the epidemics
PlotCosts <- function(simoutput,add=FALSE, ylim=NULL, ...)
  {
    if(add)
      {
        lines(simoutput$C,lty=2,lwd=2)
      }
    else
      {
        if(is.null(ylim)) ylim=c(min(simoutput$C),max(simoutput$C))
        plot(simoutput$C,type='l',lty=1,lwd=2,xlab="time",ylab="cost",
             bty='n',ylim=ylim, ...)
      }
  }
    
## Function for plotting density kernels for the sampled parameter values alongside
## the prior densities.  Arguments are the sampled values, the vectors of hyperparameters,
## and the prior density functions for each parameter.  If the
## true values of the parameters are provided, a point is placed on the x axis to indicate the
## true value on the plots. 
PlotParams <-  function(samp, dens, true, hyper=NULL) 
  {
    allnames2 <- c(expression(b),expression(k),expression(nu),expression(mu))
    allnames <- c("b","k","nu","mu")
    par(mfrow=c(2,2))
    for(p in 1:4)
      {    
        p.dens <- density(samp[,p],from=0)
        x <- p.dens$x[p.dens$x > 0]

        if(! is.null(dens)){
          p.prior.dens <- cbind(x, dens[[p]](x) )
          myylim <- c(0,max(p.dens$y, p.prior.dens[,2]))
        }  else{
          myylim <- c(0,max(p.dens$y))
        }

        if(!is.null(hyper)) {
          h <- hyper[[paste(allnames[p], "h", sep="")]]
          if(allnames[p] == "b" || allnames[p] == "k") {
            prior <- rgamma(10000, h[1], scale=h[2])
          } else {
            prior <- rbeta(10000, h[1], h[2])
            prior <- -log(1-prior)
          }
          prior.dens <- density(prior,from=0)
          myylim <- c(0,max(prior.dens$y,p.dens$y))
        }
        plot(p.dens,bty='n',xlab=allnames2[p],ylab="density",main="",ylim=myylim,type='n') 
                
        q1 <- quantile(samp[,p],.025)
        m <- mean(samp[,p])
        q3 <- quantile(samp[,p],.975)
        xvals <- c(p.dens$x[p.dens$x>=q1 & p.dens$x<=q3])
        yvals <- c(p.dens$y[p.dens$x>=q1 & p.dens$x<=q3])
        polygon(c(xvals,rev(xvals)),c(yvals,rep(0,length(xvals))),col=rgb(.8,.8,.8),border=NA)
        lines(p.dens,lwd=2)

        if(! is.null(dens)) lines(p.prior.dens,lwd=2,col=rgb(.3,.3,.3))
        if(!is.null(hyper)){
          lines(prior.dens,col=2,lwd=2)
        }
        
        points(m,0,pch='X')
        
        if(!is.null(true))
          points(true[[allnames[p]]],0,pch=19)
      }    
  }

## takes a matrix in which each row is a new time point, and the columns are sampled values of some variable
TimeSeriesOfDensities <- function(X,times,xbounds,varname="",main="")
  {
    t <- nrow(X)
    n <- ncol(X)
    plotvals <- array(0,c(t,2,512))
    for(i in 1:t)
      {
        dens <- density(X[i,],from=xbounds[1],to=xbounds[2])
        plotvals[i,2,] <- dens$x
        plotvals[i,1,] <- times[i]+ dens$y / max(dens$y) * 0.9
      }
    xrange <- c(min(times),max(times)+1)
    yrange <- c(min(plotvals[,2,]),max(plotvals[,2,]))
    plot(0,xlim=xrange,ylim=yrange,xlab='time',ylab=varname,type='n',main=main,axes=FALSE)
    for(i in 1:t)
      {
        lines(plotvals[i,1,],plotvals[i,2,],lwd=2)
        lines(c(times[i],times[i]),yrange,lty=2,lwd=1,col=rgb(.5,.5,.5))
      }
    axis(1,times)#,0:(t-1))
    axis(2)
    box()
  }

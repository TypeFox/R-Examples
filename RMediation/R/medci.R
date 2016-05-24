medci <-function(mu.x,mu.y,se.x,se.y,rho=0,alpha=.05,type="dop", plot=FALSE,plotCI=FALSE, n.mc=1e5,...)
{
    if(!is.numeric(mu.x))
        stop("Argument mu.x must be numeric!")
    if(!is.numeric(mu.y))
        stop("Argument mu.y must be numeric!")
    if(!is.numeric(se.x))
        stop("Argument se.x must be numeric!")
    if(!is.numeric(se.y))
        stop("Argument se.y must be numeric!")
    if(!is.numeric(alpha))
        stop("Argument alpha  must be numeric!")
    if(!is.numeric(rho))
        stop("Argument rho  must be numeric!")
    if(alpha<=0 || alpha >=1)
        stop("alpha must be between 0 and 1!")
    if(rho<=-1 || rho>=1)
        stop("rho must be between -1 and 1!")
    if(!is.numeric(n.mc) || is.null(n.mc))
        n.mc=1e5 # sets n.mc to default
    if (plot==TRUE)
      {
        mean.v <- c(mu.x,mu.y)
        var.mat <- matrix(c(se.x^2,se.x*se.y*rho,se.x*se.y*rho,se.y^2),2)
        x_y <- matrix(rnorm(2*n.mc),ncol=n.mc)
        x_y <- crossprod(chol(var.mat),x_y)+mean.v
        x_y <- t(x_y)
        xy <- x_y[,1]*x_y[,2]
        se.xy <- sqrt(se.y^2*mu.x^2+se.x^2*mu.y^2+2*mu.x*mu.y*rho*se.x*se.y+se.x^2*se.y^2+se.x^2*se.y^2*rho^2);
        mu.xy <- mu.x*mu.y+rho*se.x*se.y
        max1<-mu.xy+6*se.xy
        min1<-mu.xy-6*se.xy
        if (min1>0 || max1<0 )
          xrange<-round(seq(min1,max1,length=7),1)
	else
          xrange<-round(cbind(seq(min1,0,length=3),seq(0,max1,length=3)),1)
        xy <- xy[xy>min1 & xy<max1]
        plot(density(xy),xlab=expression(paste("Product ", (italic(xy)))),ylab="Density",axes=FALSE,xlim=c(min1,max1),main="",...)
        axis(1,xrange);axis(2)
        smidge <- par("cin")*abs(par("tcl"))
        text(max1-(max1-min1)/7,(par("usr")[4]),pos=1, bquote(mu== .(round(mu.xy,3)) ),...)
        text(max1-(max1-min1)/7,(par("usr")[4]-1.5*par("cxy")[2]),pos=1, bquote(sigma== .(round(se.xy,3)) ),...)
        if(plotCI){
          yci<-par("usr")[3]+diff(par("usr")[3:4])/25
          yci<-0
          MedCI <- medciMeeker(mu.x, mu.y, se.x, se.y, rho, alpha)
          arrows(MedCI[[1]][1],yci,MedCI[[1]][2],yci,length=smidge,angle=90,code=3,cex=1.5,...)
          points(mu.xy,yci,pch=19,cex=1.5,...)
          text(max1-(max1-min1)/7,(par("usr")[4]-3*par("cxy")[2]),pos=1, paste("LL=",round(MedCI[[1]][1],3)),...)
          text(max1-(max1-min1)/7,(par("usr")[4]-4.5*par("cxy")[2]),pos=1, paste("UL=",round(MedCI[[1]][2],3)),...)
        }
      }

    if (type=="all" || type=="All" || type=="ALL")
      {
        MCCI= medciMC(mu.x, mu.y, se.x, se.y, rho , alpha , n.mc = n.mc)
        asympCI <- medciAsymp(mu.x, mu.y, se.x, se.y, rho, alpha) # added 3/28/14-DT
        res <- list( MeekerCI, MCCI, asympCI)
        names(res) <- c( "Monte Carlo", "Asymptotic Normal")
        return(res)
      }
    else if (type=="DOP" || type=="dop" || type=="prodclin")
      {
        MeekerCI=medciMeeker(mu.x, mu.y, se.x, se.y, rho, alpha)
        return(MeekerCI)
      }
    else if (type=="MC" || type=="mc" || type=="Mc")
    {
      MCCI= medciMC(mu.x, mu.y, se.x, se.y, rho , alpha , n.mc = n.mc)
      return(MCCI)
    }
    else if(type=="Asymp" || type=="asymp")
      {
        asympCI <- medciAsymp(mu.x, mu.y, se.x, se.y, rho, alpha) # Modified/ added 3/28/14
        return(asympCI)
      }
    else stop("Wrong type! please specify type=\"all\", \"DOP\", \"prodclin\",\"MC\", or \"asymp\" ")
  }

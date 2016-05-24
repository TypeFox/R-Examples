############################################
#Scatter plot and box plot summary on a DiagTest3Grp object
#########################################
plot.DiagTest3Grp <- function(x,...)
  {
    ####input:
    ###(1)x: a DiagTest3Grp object
    ###(2)...: other arguments for plot function
    ###Output: scatter plot and box plot
    marker <- c(x$dat$x,x$dat$y,x$dat$z)
    grp <- factor(rep(c("D-","D0","D+"),x$dat.summary$n),levels=c("D-","D0","D+"))
    
    nobs <- length(marker)
    if(nobs!=length(grp)) stop("length(marker)!=length(grp)")
    
    par(mfrow=c(1,2))
    ##scatter plot
    plot(marker,col=as.character(factor(grp,labels=c("green","blue","red"))),...)
    abline(h=x$cut.point[1],lty=2)
    text(0.9*nobs,x$cut.point[1]+0.15,paste("t-=",round(x$cut.point[1],2),sep=""))
    abline(h=x$cut.point[2],lty=2)
    text(0.2*nobs,x$cut.point[2]+0.15,paste("t+=",round(x$cut.point[2],2),sep=""))
    legend("bottomright",legend=c(paste(x$type,"=",round(x$estimate,2),sep=""),paste((1-x$alpha)*100,"% CI=",round(x$CI[1],2),"~",round(x$CI[2],2),sep="")),bty="n")
    ##boxplot
    boxplot(marker~factor(grp),col=c("green","blue","red"))
    abline(h=x$cut.point[1],lty=2)
    abline(h=x$cut.point[2],lty=2)    
    
  }

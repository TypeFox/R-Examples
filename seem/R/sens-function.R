# function to calculate sens metrics for multiple runs

sens <- function(t.X, param, fileout, pdfout=F) {

# sensitivity analysis
# arguments:
# t.X time seq and model results
# param list of plabel, param, and factor values
# fileout prefix to store file and plots

 # for labeling graphs
 ndex=6; label.m <- c("Diff Max", "Diff Min", "Diff End", "Avg Diff", "RMS Diff", "Max Diff")

 p <- param$pval; plabel <- param$plab; pnom <- param$pnom
 np <- length(param$pval)
 # time and variables
 t <- t.X[,1]; X <- t.X[,-1]; nt <- length(t)

# sort according to order of p
 o <- order(p)
 p <- p[o];X <- X[,o]
 # finds position of nominal in the set of parameter values 
 nnom <- which(p==pnom)

 # create arrays and set to zero
 d.max <- array()
 d.min <- array()
 d.end <- array()
 avg.d <- array()
 err.m <- matrix(NA,nt,np) 
 rms.d <- array()
 max.d <- array()

 # for each parameter value calculate the metrics
 for(j in 1:np){
   d.max[j] <- max(X[,j])
   d.min[j] <- min(X[,j])
   d.end[j] <- X[nt,j]
   avg.d[j] <- mean(X[,j])
   err.m[,j] <- X[,j] - X[,nnom]
   rms.d[j] <- sqrt(mean((err.m[,j])^2))
   max.d[j] <- max(abs(err.m[,j]))
}
 y <- cbind(d.max, d.min, d.end, avg.d, rms.d, max.d)
 val.val <- data.frame(round(cbind(p,y),3))
 names(val.val) <- c("Param Val",label.m) 

 err.r <- err.m
 for(j in 1:np){
   for(i in 1:nt){
    if (X[i,nnom] ==0) err.r[i,j] <- 0 else
      err.r[i,j] <- (err.m[i,j])/X[i,nnom]
   }
   rms.d[j] <- sqrt(mean((err.r[,j])^2))
   max.d[j] <- max(abs(err.r[,j]))
 }
 yp <- cbind(d.max, d.min, d.end, avg.d, rms.d, max.d)
 
 # % change in parameter (abs is needed for negative values)
 pn <- p[nnom]; pp <- 100*(p - pn)/abs(pn)

 # % change in metrics
 yn <- y[nnom,]
  for(i in 1:4){
 # checks for all zero values for a metric
   if(yn[i] == 0) yp[,i]=0
 # calculate % if nonzero
   else yp[,i] <- 100*(y[,i]-yn[i])/yn[i] 
  }
  for(i in 5:6) {
    for(j in 1:np){
       yp[j,i] <- 100*(yp[j,i]) 
       if(j <= nnom) yp[j,i] <- -yp[j,i] 
    }
  }

 perc.perc <- data.frame(round(cbind(pp,yp),3))
 names(perc.perc) <- c("Param Perc",label.m) 

 #ratio of % changes to % parameter
 ypp <- y; ypp[,]=0
 for(i in 1:length(p))
 ypp[i,] <- yp[i,]/pp[i] 
 #interpolate nominal because div by zero
 ypp[nnom,] <-(ypp[nnom+1,]-ypp[nnom-1,])/(pp[nnom+1]-pp[nnom-1])*
              (pp[nnom]-pp[nnom-1]) + ypp[nnom-1,]

 perc.ratio <- data.frame(round(cbind(pp,ypp),3))
 names(perc.ratio) <- c("Param Perc",label.m) 

 # regression for random sampling
 if(param$fact==F){
 ycoeff <- matrix(nrow=ndex,ncol=2); yp.est <- yp;r2 <- array()
 for(i in 1:4){
  ycoeff[i,] <- lm(yp[,i]~0+pp)$coeff
  r2[i] <- summary(lm(yp[,i]~0+pp))$r.square
  yp.est[,i] <- ycoeff[i,1]*pp 
 }
 for(i in 5:6){
  #ppa <- abs(pp)
  ppa <- pp
  ycoeff[i,] <- lm(yp[,i]~0+ppa)$coeff
  r2[i] <- summary(lm(yp[,i]~0+ppa))$r.square
  yp.est[,i] <- ycoeff[i,1]*ppa 
 }

 slope <- ycoeff[,1]
 reg.slope.R2 <- data.frame(round(rbind(slope, r2),3))
 names(reg.slope.R2) <- label.m 
 } 

 # graphics
 if(pdfout==T) pdf(paste(fileout,".pdf",sep=""))

 mat <- matrix(1:2,2,1,byrow=T)
 nf <- layout(mat, widths=rep(7,2), heights=rep(7/2,2), TRUE)
 par(mar=c(4,4,1,.5),xaxs="r",yaxs="r")

 # plot val vs val
 vx <- p
 xext <- max(vx)+ 0.2*abs(max(vx)-min(vx)); xextl <- max(vx)+ 0.05*abs(max(vx)-min(vx))
 matplot(p, y, type="b", ylim=c(min(y),max(y)), xlim=c(min(p),xext),
          xlab=plabel, ylab="Metric", lty=1:ndex, pch=1:ndex,col=1,lwd=1.3)
 legend (xextl, max(y), legend=label.m, lty=c(1:ndex), pch=1:ndex, col=1,lwd=1.3,cex=0.7)
 # plot perc vs perc
 vx <- pp
 xext <- max(vx)+ 0.2*abs(max(vx)-min(vx)); xextl <- max(vx)+ 0.05*abs(max(vx)-min(vx))
 matplot(pp, yp, type="b", ylim=c(min(yp),max(yp)), xlim=c(min(pp),xext),
          xlab=paste(plabel,"%"), ylab="Metric %",  lty=1:ndex, pch=1:ndex,col=1,lwd=1.3)
 legend (xextl, max(yp), legend=label.m,  lty=1:ndex, pch=1:ndex,col=1,lwd=1.3,cex=0.7)

 if(pdfout==F){
 win.graph()
  mat <- matrix(1:2,2,1,byrow=T)
  nf <- layout(mat, widths=rep(7,2), heights=rep(7/2,2), TRUE)
  par(mar=c(4,4,1,.5),xaxs="r",yaxs="r")
 }

 # plot the perc ratio vs perc
 vx <- pp
 xext <- max(vx)+ 0.2*abs(max(vx)-min(vx)); xextl <- max(vx)+ 0.05*abs(max(vx)-min(vx))
 matplot(pp, ypp, type="b", ylim=c(min(ypp), max(ypp)), xlim=c(min(pp),xext),
         xlab=paste(plabel,"%"), ylab="Sens %/%",lty=1:ndex, pch=1:ndex,col=1,lwd=1.3)
 legend (xextl, max(ypp), legend=label.m, lty=1:ndex, pch=1:ndex, col=1,lwd=1.3,cex=0.7)
 #text(0,max(ypp),"Ratio at p=0% is indeterminate",cex=0.6) 

 # plot regression for random sampling
 if(param$fact==F){
 vx <- pp
 xext <- max(vx)+ 0.2*abs(max(vx)-min(vx)); xextl <- max(vx)+ 0.05*abs(max(vx)-min(vx))
  plot(pp,yp[,1], xlim=c(min(pp),xext),ylim=c(min(pp[1]*ycoeff[,1]),max(yp)),
       xlab=paste(plabel,"%"), ylab="Metric %", pch=1,col=1,lwd=1.3)
  lines(c(pp[1],pp[length(p)]),ycoeff[1,1]*c(pp[1],pp[length(p)]),lty=1,lwd=1.3)

 for(i in 2:ndex){
  points(pp,yp[,i],pch=i)
  lines(c(pp[1],pp[length(p)]),ycoeff[i,1]*c(pp[1],pp[length(p)]),lty=i,lwd=1.3)
 }
 legend (xextl, max(yp), legend=label.m, lty=1:ndex, pch=1:ndex, col=1,lwd=1.3,cex=0.7)
 }

 if(pdfout==T) dev.off()

 if(param$fact==T) return(list(val.val=val.val, perc.perc=perc.perc, perc.ratio=perc.ratio))
 else return(list(val.val=val.val, perc.perc=perc.perc, perc.ratio=perc.ratio, reg.slope.R2=reg.slope.R2))
}


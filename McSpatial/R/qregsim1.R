qregsim1 <- function(formall, formx, bmat, taumat, xvalues=NULL, ytarget=NULL, xcolors=NULL,graphx=TRUE, graphy=TRUE, graphsim=TRUE,  
  histogram=FALSE, histfreq=FALSE, yname=NULL, xname=NULL, nsim=0, bwadjust=1, legloc="topright",data=NULL)  {

  xlabel = xname
  ylabel = yname

  xmat <- model.frame(formall,data=data)
  y <- xmat[,1]
  yname = colnames(xmat)[1]
  if (identical(data,NULL)) {data <- xmat}

  xmat <- model.matrix(formx,data=data)
  xname <- colnames(model.frame(formx,data=data))
  x <- as.numeric(data[,xname])


  xmat <- model.matrix(formall,data=data)
  ntau = length(taumat)
  n = nrow(xmat)
  nk = ncol(xmat)

  if (identical(xlabel,NULL)) {xlabel = xname}
  if (identical(ylabel,NULL)) {ylabel = yname}

  if (graphx==TRUE) {
    if (histogram==FALSE) {
      fit <- density(x,adjust=bwadjust)
      plot(fit$x,fit$y,xlab=xlabel,ylab="Density",type="l")
    }
    if (histogram==TRUE) {
      hist(x,freq=histfreq,xlab=xlabel,main="")
    }
  }

  xobs <- sample(seq(1,n),nsim,replace=TRUE)
  bobs <- sample(seq(1,ntau),nsim,replace=TRUE)
  bdim = length(dim(bmat))

  if (bdim==2) {
    b <- bmat[,1]
    ij <- bobs
  }
  if (bdim==3) {
    b <- bmat[,,1]
    ij <- (bobs-1)*n + xobs
  }
  if (nsim>0) {
    xbhat <- b[ij]
    for (j in seq(2,nk)) {
      if (bdim==2){b <- bmat[,j]}
      if (bdim==3){b <- bmat[,,j]}
      xbhat <- xbhat + b[ij]*xmat[xobs,j]
    }
   }

  if (nsim==0&bdim==2) {xbhat <- xmat%*%t(bmat) }
  if (nsim==0&bdim==3) {
    xbhat <- b
    for (j in seq(2,nk)) {
      xbhat <- xbhat + xmat[,j]*bmat[,,j]
    }
  }


  if (identical(ytarget,NULL)) {
    ytarget <- seq(min(y),max(y),length=512)
  }
  tmin = min(ytarget)
  tmax = max(ytarget)
  nt = length(ytarget)

  fit <- density(y,adjust=bwadjust,from=tmin,to=tmax,n=nt)
  densy1 <- fit$y
  densy2 <- density(xbhat,adjust=bwadjust,from=tmin,to=tmax,n=nt)$y
  if (graphy==TRUE) {
    ymin = min(densy1,densy2)
    ymax = max(densy1,densy2)
    plot(ytarget,densy1,type="l",xlab=ylabel,ylab="Density",main="Density Functions for Actual and Predicted Values",ylim=c(ymin,ymax))
    lines(ytarget,densy2,col="red")
    legend(legloc,c("Actual","Predicted"),col=c("black","red"),lwd=1) 
  }

  if (identical(xvalues,NULL)) {
    xvalues <- quantile(x,c(.25,.75))
  }
  nx = length(xvalues)
  j = which(colnames(xmat)==xname)
  if (bdim==2){b <- bmat[,j]}
  if (bdim==3){b <- bmat[,,j]}
  if (nsim>0) {xbhat <- xbhat - b[ij]*xmat[xobs,j] }
  if (nsim==0&bdim==2) {xbhat <- xbhat - as.matrix(xmat[,j])%*%t(as.matrix(bmat[,j])) }
  if (nsim==0&bdim==3) {xbhat <- xbhat - xmat[,j]*bmat[,,j]}


  densyhat <- array(0,dim=c(nt,nx)) 
  for (jj in seq(1,nx)) {
    xmat[,j] <- xvalues[jj]
    if (nsim>0) {dhat <- xbhat + b[ij]*xmat[xobs,j] }
    if (nsim==0&bdim==2) {dhat <- xbhat + as.matrix(xmat[,j])%*%t(as.matrix(bmat[,j])) }
    if (nsim==0&bdim==3) {dhat <- xbhat + xmat[,j]*bmat[,,j]}
    densyhat[,jj] <- density(dhat,adjust=bwadjust,from=tmin,to=tmax,n=nt)$y 
  }     
   

  library(RColorBrewer)
  if (identical(xcolors,NULL)) {
    if (nx==1){xcolors = "black"}
    if (nx==2){xcolors <- c("black","red")}
    if (nx>=3){xcolors = brewer.pal(nx,"Blues")}
  }
  if (length(xcolors)==1&nx>=3) {xcolors = brewer.pal(nx,xcolors)}
  
  if (graphsim==TRUE) {
    ymin = min(densyhat)
    ymax = max(densyhat)
    lvect <- as.character(xvalues)
    plot(ytarget,densyhat[,1],type="l",xlab=ylabel,ylab="Density",ylim=c(ymin,ymax),col=xcolors[1])
    if (nx>1) {
      for (jj in seq(2,nx)) {
        lines(ytarget,densyhat[,jj],col=xcolors[jj])
      }
      legend(legloc,lvect,col=xcolors,lwd=1)
    }
  }
  
  out <- list(ytarget,densyhat,densy1,densy2)
  names(out) <- c("ytarget","densyhat","densy1","densy2")
  return(out)

}

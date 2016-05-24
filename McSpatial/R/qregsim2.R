qregsim2 <- function(formall, formx, dataframe1, dataframe2, bmat1, bmat2,  
  graphx=TRUE, graphb=TRUE, graphy=TRUE, graphdy=TRUE, nbarplot=10, yname=NULL, xnames=NULL, timenames=c("1","2"),
  leglocx="topright",leglocy="topright",leglocdy="topright",nsim=20000, bwadjx=1,bwadjy=1,bwadjdy=1)  {

  xmat1 <- model.matrix(formx,data=dataframe1)
  xmat2 <- model.matrix(formx,data=dataframe2)
  xcolnames <- colnames(xmat1)
  if (identical(xnames,NULL)) {xnames <- colnames(xmat1[,-1]) }
  xmat1 <- as.matrix(xmat1[,-1])
  xmat2 <- as.matrix(xmat2[,-1])
  colnames(xmat1) <- xcolnames[-1]
  colnames(xmat2) <- xcolnames[-1]

  allmat1 <- model.frame(formall,data=dataframe1)
  allmat2 <- model.frame(formall,data=dataframe2)
  y1 <- allmat1[,1]
  y2 <- allmat2[,1]
  if (identical(yname,NULL)) {yname <- colnames(allmat1)[1]}

  allmat1 <- model.matrix(formall,data=dataframe1)
  allmat2 <- model.matrix(formall,data=dataframe2)

  taumat <- as.numeric(rownames(bmat1))
  ntau = nrow(bmat1)
  n = nrow(xmat1)

  smallbmat1 <- bmat1[,colnames(xmat1)]
  smallbmat2 <- bmat2[,colnames(xmat2)]

  nk = ncol(xmat1)
  if (length(leglocx)==1){leglocx <- array(leglocx,dim=nk)}

  if (graphx==TRUE) {
    for (j in seq(1,nk)) {
      x1 <- xmat1[,j]
      x2 <- xmat2[,j]
      r1 <- table(x1)
      r2 <- table(x2)
      if (length(r1)>nbarplot|length(r2)>nbarplot) {
        xmin = min(x1,x2)
        xmax = max(x1,x2)
        fit1 <- density(x1,adjust=bwadjx,from=xmin,to=xmax)
        fit2 <- density(x2,adjust=bwadjx,from=xmin,to=xmax)
        h = mean(fit1$bw, fit2$bw)
        x <- seq(xmin,xmax,length=512)
        fit1 <- density(x1,bw=h,from=xmin,to=xmax)
        fit2 <- density(x2,bw=h,from=xmin,to=xmax)
        ymin = min(fit1$y,fit2$y)
        ymax = max(fit1$y,fit2$y)
        plot(x,fit1$y,xlab=xnames[j],ylab="Density",type="l",ylim=c(ymin,ymax))
        lines(x,fit2$y,col="red",lty="dashed")
        legend(leglocx[j],timenames,col=c("black","red"),lty=c("solid","dashed"),lwd=1)
      }
      if (length(r1)<=nbarplot&length(r2)<=nbarplot) { 
        rname1 <- rownames(r1)
        rname2 <- rownames(r2)
        dname1 <- setdiff(rname2,rname1)
        dname2 <- setdiff(rname1,rname2)
        n1 = length(dname1)
        n2 = length(dname2)
        if (n1>0) {
          d1 <- array(0,dim=n1)
          rownames(d1) <- dname1
          r1 <- c(r1,d1)
          o <- order(names(r1))
          r1 <- r1[o]
        }
        if (n2>0) {
          d2 <- array(0,dim=n2)
          rownames(d2) <- dname2
          r2 <- c(r2,d2)
          o <- order(names(r2))
          r2 <- r2[o]
        }
        r1 <- rbind(r1,r2)
        barplot(r1,beside=TRUE,legend.text=timenames,names.arg=colnames(r1),xlab=xnames[j],args.legend=list(x=leglocx[j]))
      }
    }
  }

  if (graphb==TRUE) {
    for (j in seq(1,nk)) {
      if (nk>1) {
        b1 <- smallbmat1[,j]
        b2 <- smallbmat2[,j]
      }
      if (nk==1) {
        b1 <- smallbmat1
        b2 <- smallbmat2
      }
      ymin = min(b1,b2)
      ymax = max(b1,b2)
      plot(taumat,b1,main=xnames[j],xlab="Quantile",ylab="coefficient",type="l",ylim=c(ymin,ymax))
      lines(taumat,b2,col="red",lty="dashed")
    }
  }

  
  n1 = nrow(xmat1)
  n2 = nrow(xmat2)
  xobs1 <- sample(seq(1:n1),nsim,replace=TRUE)
  xobs2 <- sample(seq(1:n2),nsim,replace=TRUE)
  bobs <- sample(seq(1:ntau),nsim,replace=TRUE)
  xhat1 <- allmat1[xobs1,]
  xhat2 <- allmat2[xobs2,]
  znames <- setdiff(colnames(allmat1),colnames(xmat1))
  if (identical(znames,"(Intercept)")) {xhat12 <- xhat1}
  if (!identical(znames,"(Intercept)")) {
    xhat12 <- cbind(xhat2[,znames],xhat1[,colnames(xmat1)])
    colnames(xhat12) <- c(znames,colnames(xmat1))
    xhat12 <- xhat12[,colnames(allmat1)]
  }
  bhat1 <- bmat1[bobs,]
  bhat2 <- bmat2[bobs,]

  xbhat11 <- rowSums(as.matrix(xhat1*bhat1))
  xbhat22 <- rowSums(as.matrix(xhat2*bhat2))
  xbhat12 <- rowSums(as.matrix(xhat12*bhat2))
  xmin = min(xbhat11,xbhat22,xbhat12)
  xmax = max(xbhat11,xbhat22,xbhat12)
  fit1 <- density(xbhat11,from=xmin,to=xmax,adjust=bwadjy)
  fit2 <- density(xbhat22,from=xmin,to=xmax,adjust=bwadjy)
  fit12 <- density(xbhat12,from=xmin,to=xmax,adjust=bwadjy)
  h = mean(fit1$bw,fit2$bw,fit12$bw)
  yhat11 <- density(xbhat11,from=xmin,to=xmax,bw=h)$y
  yhat12 <- density(xbhat12,from=xmin,to=xmax,bw=h)$y
  yhat22 <- density(xbhat22,from=xmin,to=xmax,bw=h)$y
  if (bwadjy==bwadjdy) {
    d2211 <- yhat22 - yhat11
    d2212 <- yhat22 - yhat12
    d1211 <- yhat12 - yhat11
  }
  if (bwadjy!=bwadjdy) {
    fit1 <- density(xbhat11,from=xmin,to=xmax,adjust=bwadjdy)
    fit2 <- density(xbhat22,from=xmin,to=xmax,adjust=bwadjdy)
    fit12 <- density(xbhat12,from=xmin,to=xmax,adjust=bwadjdy)
    h = mean(fit1$bw,fit2$bw,fit12$bw)
    temp11 <- density(xbhat11,from=xmin,to=xmax,bw=h)$y
    temp12 <- density(xbhat12,from=xmin,to=xmax,bw=h)$y
    temp22 <- density(xbhat22,from=xmin,to=xmax,bw=h)$y
    d2211 <- temp22 - temp11
    d2212 <- temp22 - temp12
    d1211 <- temp12 - temp11
  }

  ytarget <- seq(xmin,xmax,length=512)
  if (graphy==TRUE) {
    ymin = min(yhat11,yhat22,yhat12)
    ymax = max(yhat11,yhat22,yhat12)
    plot(ytarget,yhat11,type="l",xlab=yname,ylab="Density",main="Density of Quantile Predictions",ylim=c(ymin,ymax))
    lines(ytarget,yhat22,col="red")
    lines(ytarget,yhat12,col="blue")
    legend(leglocy,c(timenames,"Counterfactual"),col=c("black","red","blue"),lwd=1)
  }
  if (graphdy==TRUE) {
    ymin = min(d2211,d2212,d1211)
    ymax = max(d2211,d2212,d1211)
    plot(ytarget,d2211,type="l",xlab=yname,ylab="Change in Density",main="Density Change Decomposition",ylim=c(ymin,ymax))
    lines(ytarget,d2212,col="red")
    lines(ytarget,d1211,col="blue")
    legend(leglocdy,c("Total","Variables","Coefficients"),col=c("black","red","blue"),lwd=1)
  }
  
  out <- list(ytarget,yhat11,yhat22,yhat12,d2211,d2212,d1211)
  names(out) <- c("ytarget","yhat11","yhat22","yhat12","d2211","d2212","d1211")
  return(out)

}


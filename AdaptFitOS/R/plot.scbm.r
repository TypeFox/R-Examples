plot.scbm=function(x, select=NULL, bands=TRUE,
         grid=50, pages=0, plot=TRUE, ylim=NULL, xlab=NULL, ylab=NULL,
         scb.lwd=1, scb.lty="dotted", shade=FALSE, shade.col=grey(0.85),
         residuals=FALSE, residuals.col="steelblue", bayes=FALSE, rug=TRUE,...){
    fit.mat <- x$fit
    if (is.null(select)) {if (is.null(x$select)) nsmooth=1:length( x$aspobject$info$pen$degree) else nsmooth=x$select} else nsmooth=select
    drv=x$drv
    m=sum(sapply(fit.mat$smooth, FUN=function(z) return(z$m[1])))
    data <- fit.mat$model
    y <- data[,1]
    n <- length(y)
    sigma2=x$aspobject$fit$sigma^2
    Pen<-sigma2/diag(x$aspobject$aux$random.var)

    if (is.null(ylim)) ylim=vector("list",max(nsmooth))
    if (is.null(xlab)) xlab=vector("list",max(nsmooth))
    if (is.null(ylab)) ylab=vector("list",max(nsmooth))
    if (!is.list(ylim)){ylim.all=ylim; ylim=vector("list",max(nsmooth)); for (j in nsmooth) ylim[[j]]=ylim.all}
    if (!is.list(xlab)){xlab.all=xlab; xlab=vector("list",max(nsmooth)); for (j in nsmooth) xlab[[j]]=xlab.all}
    if (!is.list(ylab)){ylab.all=ylab; ylab=vector("list",max(nsmooth)); for (j in nsmooth) ylab[[j]]=ylab.all}

    plotresiduals=residuals
    if (drv>0&residuals==T) plotresiduals=F
  
    #layout as in package mgcv by Simon Wood via pages
      # sort out number of pages and plots per page
      if (pages>length(nsmooth)) pages=length(nsmooth)
      if (pages<0) pages=0
      if (pages!=0){
        ppp=length(nsmooth)%/%pages
        if (length(nsmooth)%%pages!=0){
          ppp<-ppp+1
          while (ppp*(pages-1)>=length(nsmooth)) pages<-pages-1
        }
        # now figure out number of rows and columns
        c<-trunc(sqrt(ppp))
        if (c<1) c<-1
        r<-ppp%/%c
        if (r<1) r<-1
        while (r*c<ppp) r<-r+1
        while (r*c-ppp >c && r>1) r<-r-1
        while (r*c-ppp >r && c>1) c<-c-1
        oldpar<-par(mfrow=c(r,c))
      }
      else { ppp<-1;oldpar<-par()}

      if ((pages==0&&prod(par("mfcol"))<length(nsmooth)&&dev.interactive()) || pages>1&&dev.interactive()) ask <- TRUE else ask <- FALSE
      if (ask){
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
      }

    ucb=lcb=fit=grid.x=Xj=Stdev.fit=residuals=list()
    for (j in nsmooth){
      xx <- data[,1+j]
      grid.x[[j]] <- seq(min(xx),max(xx),length=grid)
      data.grid <- data.frame(x=grid.x[[j]])
      names(data.grid) <- names(fit.mat$model)[1+j]


      if (drv==0){
        CZj.grid= Predict.matrix.lme(fit.mat$smooth[[j]],data.grid,center=F)
        if (class(fit.mat$smooth[[j]])=="ospline.smooth") Cj.grid= CZj.grid$C[,drop=F]  else Cj.grid= CZj.grid$C[,-1,drop=F]
        Zj.grid= CZj.grid$Z
        if (fit.mat$smooth[[j]]$m[1]==1) Cj.grid=matrix(0,nrow(Zj.grid),0)
        Xj.grid=  cbind(Cj.grid,Zj.grid)

        # centering with respect to curve at data points
          data1 <- data.frame(x= data[,1+j])
          names(data1) <- names(fit.mat$model)[1+j]
          CZj1= Predict.matrix.lme(fit.mat$smooth[[j]],data1,center=F)
          if (class(fit.mat$smooth[[j]])=="ospline.smooth") Cj1= CZj1$C[,drop=F]  else Cj1= CZj1$C[,-1,drop=F]
          Zj1= CZj1$Z
          if (fit.mat$smooth[[j]]$m[1]==1) Cj1=matrix(0,nrow(Zj1),0)
          Xj1=  cbind(Cj1,Zj1)

          Xj.grid=t(apply(Xj.grid,1,function(z) z-colSums(Xj1)/n))
      }
      else {
        if (class(fit.mat$smooth[[j]])=="ospline.smooth") {
          CZj.grid <- Predict.matrix.lme(fit.mat$smooth[[j]],data.grid,drv=drv,center=F)
          Cj.grid= CZj.grid$C[,drop=F]
          Zj.grid= CZj.grid$Z
          Xj.grid=  cbind(Cj.grid,Zj.grid)
        }
        else  Xj.grid=drvbasis(grid.x[[j]],degree=fit.mat$smooth[[j]]$m[1],knots=fit.mat$smooth[[j]]$knots,drv=drv,basis=class(fit.mat$smooth[[j]]))[,-1]
      }

#     Sj= Xj.grid%*%x$cov.coef[[j]]%*%t(Xj)%*%t(x$W[[j]])
      Sj= Xj.grid%*%x$cov.coef[[j]]%*%t(x$WjXj[[j]])
      
      if (bayes==F){
          if (!x$hetero) SS <- tcrossprod(Sj)
          else SS <- (Sj)%*%diag(x$sigmax.fitted^2)%*%t(Sj)
          Stdev.fit[[j]] <- sqrt(diag(SS)*(sigma2))
      }  
      else{
          S= Xj.grid%*%x$cov.coef[[j]]%*%t(Xj.grid)
          Stdev.fit[[j]] <- sqrt(diag(S)*(sigma2))
      }
      fit[[j]]=Sj%*%y
      ucb[[j]] <- fit[[j]]+x$crit[[j]]*Stdev.fit[[j]]
      lcb[[j]] <- fit[[j]]-x$crit[[j]]*Stdev.fit[[j]]

      if (plot) {
        if (plotresiduals)  residuals[[j]] <-  x$fitted[[j]]+ (y-x$aspobject$fitted)
        if (is.null(ylim[[j]])){
          if (!plotresiduals) ylim[[j]]=c(min(lcb[[j]]),max(ucb[[j]]))
          else         ylim[[j]]=c(min(residuals[[j]]),max(residuals[[j]]))
        }
        if (is.null(ylab[[j]])) ylab[[j]]=paste("s", paste(rep("\'",drv),collapse="","",sep=""),"(",names(fit.mat$model)[j+1],")",sep="")
        if (is.null(xlab[[j]])) xlab[[j]]=names(fit.mat$model)[j+1]

        plot(grid.x[[j]],fit[[j]],type="l",ylim= ylim[[j]],xlab=xlab[[j]],ylab=ylab[[j]],...)#

        if (shade){
           polygon(c(grid.x[[j]], rev(grid.x[[j]])), c(lcb[[j]], rev(ucb[[j]])), col = shade.col, border = NA)
           lines(grid.x[[j]],fit[[j]],...)
        }
        if (plotresiduals){points(data[,1+j], residuals[[j]], col=residuals.col); lines(grid.x[[j]],fit[[j]],...)}

        if (rug) rug(xx)
        if(bands==TRUE){
            lines(grid.x[[j]],lcb[[j]],lty=scb.lty,lwd=scb.lwd)
            lines(grid.x[[j]],ucb[[j]],lty=scb.lty,lwd=scb.lwd)
        }
        if (drv>0) abline(h=0)
       }
     }
 names(grid.x)=names(fit)=names(lcb)=names(ucb)=names(Stdev.fit)=(names(fit.mat$model)[-1])[nsmooth]

if (pages>0) par(oldpar)
invisible(list(grid.x=grid.x,fitted=fit,lcb=lcb,ucb=ucb,drv=drv,Stdev.fit=Stdev.fit,ylim=ylim,residuals=residuals))
}

`plot.scbTest` <-
function(x,bands=TRUE,grid=50, plot=TRUE,mfrow=NULL,select=NULL,ylim=NULL,xlab=NULL,ylab=NULL,scb.lwd=1,scb.lty="dotted",shade=FALSE, shade.col=grey(0.85),residuals=FALSE,residuals.col="steelblue",rug=TRUE,...)   # hetero=F,
  {
    fit.mat <- x$fit
    if (is.null(select)) {if (is.null(x$select)) nsmooth=1:length( x$aspobject$info$pen$degree) else nsmooth=x$select} else nsmooth=select
    drv=x$drv
    m=sum(sapply(fit.mat$smooth, FUN=function(z) return(z$m[1])))
    data <- fit.mat$model
    y <- data[,1]
    n <- length(y)
    sigma2=x$aspobject$fit$sigma^2

    if (is.null(ylim)) ylim=vector("list",max(nsmooth))
    if (is.null(xlab)) xlab=vector("list",max(nsmooth))
    if (is.null(ylab)) ylab=vector("list",max(nsmooth))
    if (plot) if(all(par()$mfrow==c(1,1))){if (is.null(mfrow)) par(mfrow=c(1,length(nsmooth))) else par(mfrow=mfrow)}

    ucb=lcb=fit=grid.x=Xj=Stdev.fit=residuals=list()

    for (j in nsmooth){
      xx <- data[,1+j]
      grid.x[[j]] <- seq(min(xx),max(xx),length=grid)
      data.grid <- data.frame(x=grid.x[[j]])
      names(data.grid) <- names(fit.mat$model)[1+j]


    if (drv==0){
      CZj.grid= Predict.matrix.lme(fit.mat$smooth[[j]],data.grid)
      if (class(fit.mat$smooth[[j]])=="ospline.smooth") Cj.grid= CZj.grid$C[,drop=F]  else Cj.grid= CZj.grid$C[,-1,drop=F]
      Zj.grid= CZj.grid$Z
   #   Xj.grid=  cbind(Cj.grid,Zj.grid)
    }
    else {
      if (class(fit.mat$smooth[[j]])=="ospline.smooth") {
        CZj.grid <- Predict.matrix.lme(fit.mat$smooth[[j]],data.grid,drv=drv,center=F)
        Cj.grid= CZj.grid$C[,drop=F]
        Zj.grid= CZj.grid$Z
    #    Xj.grid=  cbind(Cj.grid,Zj.grid)
      }
      else  Xj.grid=drvbasis(grid.x[[j]],degree=fit.mat$smooth[[j]]$m[1],knots=fit.mat$smooth[[j]]$knots,drv=drv,basis=class(fit.mat$smooth[[j]]))[,-1]
    }



#     Sj= Xj.grid%*%x$cov.coef[[j]]%*%t(Xj)%*%t(x$W[[j]])
      Sj= Zj.grid%*%x$cov.coef[[j]]%*%t(x$WjXj[[j]])
      if (!x$hetero) SS <- tcrossprod(Sj)
      else SS <- (Sj)%*%diag(x$sigmax.fitted^2)%*%t(Sj)

      fit[[j]]=Sj%*%y
      Stdev.fit[[j]] <- sqrt(diag(SS)*(sigma2))
      ucb[[j]] <- fit[[j]]+x$crit[[j]]*Stdev.fit[[j]]
      lcb[[j]] <- fit[[j]]-x$crit[[j]]*Stdev.fit[[j]]

      if (plot) {

      #   if (residuals==T)  residuals[[j]] <-  x$fitted[[j]]+ (y-x$aspobject$fitted)

        if (is.null(ylim[[j]])){
          ylim[[j]]=c(min(lcb[[j]]),max(ucb[[j]]))
        }
        if (is.null(ylab[[j]])) ylab[[j]]=paste("s", paste(rep("\'",drv),collapse="","",sep=""),"(",names(fit.mat$model)[j+1],")",sep="")
        if (is.null(xlab[[j]])) xlab[[j]]=names(fit.mat$model)[j+1]
        plot(grid.x[[j]],fit[[j]],type="l",ylim= ylim[[j]],xlab=xlab[[j]],ylab=ylab[[j]],...)#
         if (shade){
           polygon(c(grid.x[[j]], rev(grid.x[[j]])), c(lcb[[j]], rev(ucb[[j]])), col = shade.col, border = NA)
           lines(grid.x[[j]],fit[[j]],...)
        }
      #  if (residuals){points(data[,1+j], residuals[[j]], col=residuals.col); lines(grid.x[[j]],fit[[j]],...)}

        if (rug) rug(xx)
        if(bands==TRUE){
            lines(grid.x[[j]],lcb[[j]],lty=scb.lty,lwd=scb.lwd)
            lines(grid.x[[j]],ucb[[j]],lty=scb.lty,lwd=scb.lwd)
        }
        if (drv>0) abline(h=0)
        abline(h=0)
       }
     }
 names(grid.x)=names(fit)=names(lcb)=names(ucb)=names(Stdev.fit)=(names(fit.mat$model)[-1])[nsmooth]

invisible(list(grid.x=grid.x,fitted=fit,lcb=lcb,ucb=ucb,drv=drv,Stdev.fit=Stdev.fit,ylim=ylim,residuals=residuals))
}

rsplot.gls <- function (model, group="none", id=FALSE){  

    mf <- model.frame(model)
    yname <- names(attr(terms(model),"dataClasses"))[1]
    fittedv <- fitted(model)
    residv <- resid(model, type="pearson")
    obsv <- eval(parse(text=yname), mf)
    if (length(obsv)!=length(fittedv)) obsv <- na.omit(obsv)
    
    op <- par(mfrow=c(2, 2), cex=0.6, mar=c(5, 5, 4, 2), mex=0.8)
    
   	form <- attr(model$modelStruct$corStruct, "formula")
    if (is.null(form)) form <- ~1
    plot(ACF ~ lag, type="h", col="blue", data=nlme::ACF(model, form=form))
    abline(h=0)
    abline(h=0.1, col="blue", lty=2)
    abline(h=-0.1, col="blue", lty=2)
    
    tmp <- qqnorm(residv, col="blue", main = "Normal Plot for Residuals", 
      xlab="Standardized residuals",  ylab="Quantiles of Standard Normal")	
    qqline(residv)
  	lab.df <- data.frame(x=tmp$x, y=tmp$y, lab = rownames(mf))
  	lab.df <- lab.df[order(lab.df$x),]
  	with(lab.df[c(1:3, (dim(lab.df)[1]-2):dim(lab.df)[1]),], text(x,y, labels = lab))
    
    if (unique(group=="none")) {
      plot(residv~fittedv, col="blue", main = "Residuals vs Fitted", 
        xlab="Fitted values", ylab="Standardized residuals")
      loess.fit <- loess.smooth(fittedv,residv)
      lines(loess.fit$x, loess.fit$y, col="red", lty=1, lwd=0.4)
    }else{
      if(unique(group %in% names(mf))){
        colf <- as.numeric(mf[, group])
      }else{
        colf <- as.numeric(group)
        if (!length(colf)==nrow(mf)) stop("The length of 'group' must be the same as the other predicted variable!")
      }
      plot(residv~fittedv, col=colf, main = "Residuals vs Fitted", 
        xlab="Fitted values", ylab="Standardized residuals")
      loess.fit <- loess.smooth(fittedv,residv)
      lines(loess.fit$x, loess.fit$y, col="red", lty=1, lwd=0.4)
    }
    abline(0,0)
    if (id) identify(fittedv, residv, labels = rownames(mf))
    
    plot(fittedv~obsv, xlab=yname, ylab="Fitted Values", col="blue", main = "Fitted vs Observed")
    abline(0,1)

#    CookD <- CookD(model, plot=FALSE)
#	outD <- CookD >= sort(CookD, decreasing =TRUE)[3]                        # Outlying Di's
#    labid <- names(CookD)
#    plot(CookD,  xlab="Obs. number", col="blue", ylim=c(0, max(CookD)+0.005),
#         main="Cook's Distance", ylab = "Cook's distance", type = "h")
#    text((1:length(CookD))[outD], CookD[outD], labid[outD], pos=3)          # Annotation
#    points((1:length(CookD))[outD], CookD[outD], pch=16, col="blue")
    par(op)
  }

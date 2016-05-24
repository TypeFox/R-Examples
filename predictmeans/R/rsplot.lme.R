rsplot.lme <- function (model, group="none", level=1, slope=FALSE, id=FALSE){  #model@call$family

  mf <- model.frame(model)
#  tm <- terms(model) 
#  if (class(model)[1]%in%c("lmerMod", "glmerMod")) cls <- sapply(all.vars(formula(model, fixed.only=TRUE)),function(x) class(model.frame(model)[[x]])[1]) else cls <- attr(tm,"dataClasses")
#  yname <- names(cls)[1]
  yname <- as.character(attr(terms(model), "variables"))[[2]]

  if (class(model)[1]=="lme") {
    rand.str <- nlme::ranef(model)
    if (!is.data.frame(rand.str)) {
      if (slope) {
        qqy <- nlme::ranef(model)[[level]][,2]
        mtitle <- "Normal Plot for Random Slope"
      }else{
        qqy <- nlme::ranef(model)[[level]][,1]
        mtitle <- "Normal Plot for Random Intercept"
      } # end of slope if
    }else{
      if (slope) {
        qqy <- nlme::ranef(model)[,2]
        mtitle <- "Normal Plot for Random Slope"
      }else{
        qqy <- nlme::ranef(model)[,1]
        mtitle <- "Normal Plot for Random Intercept"
      } # end of slope if
    }# end of rand,str if
    fittedv <- fitted(model)
    residv <- resid(model, type="pearson")
    obsv <- eval(parse(text=yname), mf)
    if (length(obsv)!=length(fittedv)) obsv <- na.omit(obsv)
  }else{
    if (class(model)[1]%in%c("lmerMod", "glmerMod")){
      if (slope) {
        qqy <- as.data.frame(lme4::ranef(model)[[level]])[,2]
        mtitle <- "Normal Plot for Random Slope"
      }else{
        qqy <- as.data.frame(lme4::ranef(model)[[level]])[,1]
        mtitle <- "Normal Plot for Random Intercept"
      } # end of slope if
      fittedv <- na.omit(fitted(model))
      residv <- na.omit(resid(model))
      obsv <- mf[,yname]
      if (!is.null(dim(obsv))) obsv <- obsv[,1]/rowSums(obsv)
    }else{
      stop("This function only works for \'lme\' or \'lmer\' models!")
    }
  } # end of mcls if
  op <- par(mfrow=c(2, 2), cex=0.6, mar=c(5, 5, 4, 2), mex=0.8)
  qqnorm(qqy, col="blue", main = mtitle)
  qqline(qqy)
  
  if(class(model)[1]=="glmerMod") {
  ##plot random effects against the predicted values and check for no trend:
    mu <- model.matrix(model)%*%lme4::fixef(model)
    RandomEffects <- t(as.matrix(model@pp$Zt)) %*% unlist(lme4::ranef(model))
    plot(mu, RandomEffects, col = "blue", main = "Random efects vs Fitted (mu)")
    abline(0,0)
  }else{
    tmp <- qqnorm(residv, col="blue", main = "Normal Plot for Residuals",
	  xlab="Standardized residuals",  ylab="Quantiles of Standard Normal")
    qqline(residv)
	lab.df <- data.frame(x=tmp$x, y=tmp$y, lab = rownames(mf))
	lab.df <- lab.df[order(lab.df$x),]
	with(lab.df[c(1:3, (dim(lab.df)[1]-2):dim(lab.df)[1]),], text(x,y, labels = lab))
  }
  
  if (unique(group=="none")) {
    plot(residv~fittedv, col="blue", main = "Residuals vs Fitted",
	  ylab="Standardized residuals",  xlab="Fitted values")
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
	  ylab="Standardized residuals",  xlab="Fitted values")
    loess.fit <- loess.smooth(fittedv,residv)
    lines(loess.fit$x, loess.fit$y, col="red", lty=1, lwd=0.4)
  }
  abline(0,0)
  if (id) identify(fittedv, residv, labels = rownames(mf))
 
if(class(model)[1]=="glmerMod" && length(unique(model@resp$y))==2) {
## plot logistic curve, mean change, T/F 
    observed <- model@resp$y
    sf <- sort(fittedv, index=T)   # sort the fitted values
    plot(sf$x, ylim=c(0,1), type="s", col="blue", lwd=2,
    xlab="sorted sample number", ylab="probability")
    text(0, min(fittedv)+0.06,
    "fitted probability", col="blue", pos=4)
    title("Goodness of fit for logistic model")
    abline(h=mean(fittedv), lty=2)
    text(0, mean(fittedv)+.02, "mean probability", pos=4)    
    observed <- ifelse(observed==1, observed+0.025, observed-0.025)
    abline(v=length(observed)/2,lty=2)
    text(length(observed)/2,.03,"midpoint",pos=4)
    # show T/F
    points(1:length(observed), observed[sf$ix],
    pch="|",cex=1,col=ifelse(as.integer(observed[sf$ix]), "green4", "maroon1"))  
  }else{
    plot(fittedv~obsv, xlab=yname, ylab="Fitted Values", col="blue", main = "Fitted vs Observed")
    abline(0,1)
  }
  par(op)
}

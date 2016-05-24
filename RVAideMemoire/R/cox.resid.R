cox.resid <-
function (model) {
  colonnes <- colnames(model.frame(model))
  m.frame <- as.data.frame(model.frame(model)[,-1])
  names(m.frame) <- colonnes[-1]
  if (!is.null(model$xlevels)) {
    variables <- colnames(m.frame)[which(colnames(m.frame)%in%names(model$xlevels))]
    if (length(variables)<ncol(m.frame)) {
	covar <- as.data.frame(m.frame[,!colnames(m.frame)%in%variables])
	names(covar) <- colnames(m.frame)[!colnames(m.frame)%in%variables]
    } else {
	stop("no covariate in the model")
    }
  } else {
    covar <- m.frame
  }
  res <- residuals(model,type="martingale")
  if (ncol(covar)>1) {
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))
    par(mfrow=n2mfrow(ncol(covar)))
  }
  for (i in 1:ncol(covar)) {
    plot(covar[,i],res,xlab=colnames(covar)[i],ylab="Martingale residuals")
    abline(h=0,lty=3,col="grey")
    panel.smooth(covar[,i],res)
  }
}


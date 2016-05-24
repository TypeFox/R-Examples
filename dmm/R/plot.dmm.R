plot.dmm <-
function(x, traitset="all",gls=F, ...)
#  plot.dmm() - plot residuals of a dmm fitted model object
{
  if(traitset[1] == "all"){
    traits <- dimnames(x$b)[[2]][1:ncol(x$b)]
  }
  else {
    traits <- traitset
  }
  traitpairs <- permpaste(traits)

  n <- x$totn
  l <- length(traits)

  if(x$dmeopt == "qr"){   # "qr" case is different
    if(!gls) {
      if(is.null(x$dme.fit)){
        stop("Object ",substitute(x)," must contain attribute dme.fit:\n")
      }
      if(is.null(x$dme.psi)) {
        stop("Object ",substitute(x)," must contain attribute dme.psi:\n")
      }
      ymat <- matrix(x$dme.psi,n^2,l^2,dimnames=list(NULL,traitpairs))
    # histograms of residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        hist(qr.resid(x$dme.fit,x$dme.psi)[,traitpairs[i]],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Histograms of dyadic residuals for fit object ",substitute(x)))
    # qqnorm of residuals
      devAskNewPage(ask=TRUE)
      for(i in 1 : length(traitpairs)) {
        qqnorm(qr.resid(x$dme.fit,x$dme.psi)[,traitpairs[i]],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Qqnorm plots of dyadic residuals for fit object ",substitute(x)))
    # fitted values against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(qr.fitted(x$dme.fit,x$dme.psi)[,traitpairs[i]],qr.resid(x$dme.fit,x$dme.psi)[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyadic model fitted values against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],qr.resid(x$dme.fit,x$dme.psi)[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against fitted values
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],qr.fitted(x$dme.fit,x$dme.psi)[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against fitted values for fit object ",substitute(x)))
    }


    if(gls) {
      if(is.null(x$gls)) {
        stop("Object ",substitute(x)," must contain attribute gls:\n")
      }
      if(is.null(x$gls$dme.fit)){
        stop("Object ",substitute(x)," must contain attribute gls$dme.fit:\n")
      }
      if(is.null(x$gls$dme.psi)){
        stop("Object ",substitute(x)," must contain attribute gls$dme.psi:\n")
      }
      residmat <- matrix(qr.resid(x$gls$dme.fit,x$gls$dme.psi),n^2*l^2,l^2,dimnames=list(NULL,traitpairs))
      fittedmat <- matrix(qr.fitted(x$gls$dme.fit,x$gls$dme.psi),n^2*l^2,l^2,dimnames=list(NULL,traitpairs))
      ymat <- matrix(x$gls$dme.psi,n^2*l^2,l^2,dimnames=list(NULL,traitpairs))
    # histograms of residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        hist(residmat[,traitpairs[i]],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Histograms of dyadic residuals for fit object ",substitute(x)))
    # qqnorm of residuals
      devAskNewPage(ask=TRUE)
      for(i in 1 : length(traitpairs)) {
        qqnorm(residmat[,traitpairs[i]],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Qqnorm plots of dyadic residuals for fit object ",substitute(x)))
    # fitted values against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(fittedmat[,traitpairs[i]],residmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyadic model fitted values against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],residmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against fitted values
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],fittedmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against fitted values for fit object ",substitute(x)))
    }
  }

  else if(x$dmeopt == "pcr") {   # "pcr" case
    nc <- x$dme.fit$ncomp
    if(!gls) {
      if(is.null(x$dme.fit)){
        stop("Object ",substitute(x)," must contain attribute dme.fit:\n")
      }
      residmat <- array(resid(x$dme.fit),dim=c(n^2*l^2,l^2,nc),dimnames=list(NULL,traitpairs,NULL))
      fittedmat <- array(fitted(x$dme.fit),dim=c(n^2*l^2,l^2,nc),dimnames=list(NULL,traitpairs,NULL))
      ymat <- matrix(x$dme.fit$y,n^2 * l^2,l^2,dimnames=list(NULL,traitpairs))
    # histograms of residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        hist(residmat[,traitpairs[i],nc],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Histograms of dyadic residuals for fit object ",substitute(x)))
    # qqnorm of residuals
      devAskNewPage(ask=TRUE)
      for(i in 1 : length(traitpairs)) {
        qqnorm(residmat[,traitpairs[i],nc],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Qqnorm plots of dyadic residuals for fit object ",substitute(x)))
    # fitted values against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(fittedmat[,traitpairs[i],nc],residmat[,traitpairs[i],nc],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyadic model fitted values against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],residmat[,traitpairs[i],nc],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against fitted values
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],fittedmat[,traitpairs[i],nc],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against fitted values for fit object ",substitute(x)))
    }


    if(gls) {
      if(is.null(x$gls)) {
        stop("Object ",substitute(x)," must contain attribute gls:\n")
      }
      if(is.null(x$gls$dme.fit)){
        stop("Object ",substitute(x)," must contain attribute gls$dme.fit:\n")
      }
      residmat <- array(resid(x$gls$dme.fit),dim=c(n^2*l^2,l^2,nc),dimnames=list(NULL,traitpairs,NULL))
      fittedmat <- array(fitted(x$gls$dme.fit),dim=c(n^2*l^2,l^2,nc),dimnames=list(NULL,traitpairs,NULL))
      ymat <- matrix(x$gls$dme.fit$y,n^2*l^2,l^2,dimnames=list(NULL,traitpairs))
    # histograms of residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        hist(residmat[,traitpairs[i],nc],xlab=traitpairs[i],main="")
#       hist(ymat[,i],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Histograms of dyadic residuals for fit object ",substitute(x)))
    # qqnorm of residuals
      devAskNewPage(ask=TRUE)
      for(i in 1 : length(traitpairs)) {
        qqnorm(residmat[,traitpairs[i],nc],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Qqnorm plots of dyadic residuals for fit object ",substitute(x)))
    # fitted values against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(fittedmat[,traitpairs[i],nc],residmat[,traitpairs[i],nc],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyadic model fitted values against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],residmat[,traitpairs[i],nc],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against fitted values
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],fittedmat[,traitpairs[i],nc],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against fitted values for fit object ",substitute(x)))
    }
  }
  
  else {   # other than "qr" or "pcr" cases
    if(!gls) {
      if(is.null(x$dme.fit)){
        stop("Object ",substitute(x)," must contain attribute dme.fit:\n")
      }
      residmat <- matrix(resid(x$dme.fit),n^2,l^2,dimnames=list(NULL,traitpairs))
      fittedmat <- matrix(fitted(x$dme.fit),n^2,l^2,dimnames=list(NULL,traitpairs))
      ymat <- matrix(x$dme.fit$y,n^2,l^2,dimnames=list(NULL,traitpairs))
    # histograms of residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        hist(residmat[,traitpairs[i]],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Histograms of dyadic residuals for fit object ",substitute(x)))
    # qqnorm of residuals
      devAskNewPage(ask=TRUE)
      for(i in 1 : length(traitpairs)) {
        qqnorm(residmat[,traitpairs[i]],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Qqnorm plots of dyadic residuals for fit object ",substitute(x)))
    # fitted values against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(fittedmat[,traitpairs[i]],residmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyadic model fitted values against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],residmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against fitted values
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],fittedmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against fitted values for fit object ",substitute(x)))
    }


    if(gls) {
      if(is.null(x$gls)) {
        stop("Object ",substitute(x)," must contain attribute gls:\n")
      }
      if(is.null(x$gls$dme.fit)){
        stop("Object ",substitute(x)," must contain attribute gls$dme.fit:\n")
      }
      residmat <- matrix(resid(x$gls$dme.fit),n^2*l^2,l^2,dimnames=list(NULL,traitpairs))
      fittedmat <- matrix(fitted(x$gls$dme.fit),n^2*l^2,l^2,dimnames=list(NULL,traitpairs))
      ymat <- matrix(x$gls$dme.fit$y,n^2*l^2,l^2,dimnames=list(NULL,traitpairs))
    # histograms of residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        hist(residmat[,traitpairs[i]],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Histograms of dyadic residuals for fit object ",substitute(x)))
    # qqnorm of residuals
      devAskNewPage(ask=TRUE)
      for(i in 1 : length(traitpairs)) {
        qqnorm(residmat[,traitpairs[i]],xlab=traitpairs[i],main="")
      }
      mtext(side=3,line=0,cex=1.2,outer=TRUE,paste("Qqnorm plots of dyadic residuals for fit object ",substitute(x)))
    # fitted values against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(fittedmat[,traitpairs[i]],residmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyadic model fitted values against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against residuals
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],residmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against dyadic residuals for fit object ",substitute(x)))
    # observed dyad covariances against fitted values
      devAskNewPage(ask=TRUE)
      par(mfrow=c(length(traits),length(traits)),oma=c(0,0,4,0))
      for(i in 1 : length(traitpairs)) {
        plot(ymat[,traitpairs[i]],fittedmat[,traitpairs[i]],main=traitpairs[i])
      }
      mtext(side=3,line=0,cex=1,outer=TRUE,paste("Plot of dyad covariances against fitted values for fit object ",substitute(x)))
    }
  }
}

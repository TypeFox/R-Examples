"plot.MCMCglmm"<-function(x, random=FALSE, ...){

  nF<-x$Fixed$nfl
  devAskNewPage.orig<-devAskNewPage()
  if(random){
    nF<-sum(rep(x$Random$nrl, x$Random$nfl))+nF
    if(nF!=dim(x$Sol)[2]){stop("random effects not saved and cannot be plotted")}    
  }

  plot(x$Sol[,1:nF, drop=FALSE], ...)
  devAskNewPage(TRUE)
  if(is.null(x$Lambda)==FALSE){
    plot(x$Lambda, ...)
    devAskNewPage(TRUE)
  }
  plot(x$VCV, ...)
  devAskNewPage(devAskNewPage.orig)
}


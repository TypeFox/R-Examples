plot.AUCRF <-
function(x, which=c("auc","ranking","psel") ,showOpt=TRUE, digits=4, maxvars=NULL, ...){
	cl <- match.call()

  which <- match.arg(which)

  opar <- par("cex", "pch")
  on.exit(par(opar))
  n <- ifelse(is.null(maxvars),x$Kopt, maxvars)
  par(cex = max(1 - 0.1 * n %/% 10, 0.5), pch = 19)

  switch(which,
	psel={
		if(is.null(x$Psel))
			cat("No Psel information available. See AUCRFcv help.\n")
    else{
      imp <- x$Psel[x$Xopt]
	    imp <- sort(imp,decreasing=T)
	    if(!is.null(maxvars)) imp <- sort(x$Psel,decreasing=T)[1:maxvars]
      if(is.null(cl$pch)) dotchart(imp[length(imp):1],pch=par("pch"),...)
      else dotchart(imp[length(imp):1],...)
      if(is.null(cl$xlab)) title(xlab="Probability of selection")
	  }
	},
	ranking={
		imp <- x$ranking[x$Xopt]
    if(!is.null(maxvars)) imp <- sort(x$ranking,decreasing=T)[1:maxvars]
    if(is.null(cl$pch)) dotchart(imp[length(imp):1],pch=par("pch"),...)
    else dotchart(imp[length(imp):1],...)
    if(is.null(cl$xlab)) title(xlab=x$ImpMeasure)
	},
	auc={	
		par(opar)
    type <- ifelse(is.null(cl$type), "o", eval(cl$type)) 
    pch <- ifelse(is.null(cl$pch), 20, eval(cl$pch))
    col <- ifelse(is.null(cl$col), 2, eval(cl$col))

    if(is.null(cl$ylim)){
      maxAUC <- max(x$AUCcurve$AUC)
      minAUC <- min(x$AUCcurve$AUC)
      r <- maxAUC-minAUC
      ylim <- c(max(0,minAUC-r), min(1,maxAUC+r))
    }
    else{
      ylim <- eval(cl$ylim)
    }
  
    ylab <- ifelse(is.null(cl$ylab), "OOB-AUC", eval(cl$ylab))
    xlab <- ifelse(is.null(cl$xlab), "Number of selected variables", eval(cl$xlab))
     
    AUCcurve <- x$AUCcurve
    
    m <- match(c("x","wich","showOpt","digits","maxvars","type", "pch", "col", "ylim", "ylab", "xlab"), names(cl), 0L)
    plotCall <- cl[-m]
    plotCall[[1]] <- as.name("plot")
    plotCall$x <- as.name("AUCcurve")
    plotCall$type <- as.name("type")
    plotCall$pch <- as.name("pch")
    plotCall$col <- as.name("col")
    plotCall$ylim <- as.name("ylim")
    plotCall$ylab <- as.name("ylab")
    plotCall$xlab <- as.name("xlab")
    eval(plotCall)
   
    if(showOpt){
      text(x$Kopt,x$"OOB-AUCopt"+(ylim[2]-ylim[1])/10, paste("OOB-AUCopt = ",round(x$"OOB-AUCopt",digits)," (Kopt = ",x$Kopt,")",sep=""), pos=4, cex=0.7, offset=0)
      text(x$Kopt,x$"OOB-AUCopt","|", pos=3)
      if(!is.null(x$cvAUC)) 
        text(x$Kopt,x$"OOB-AUCopt"+1.5*(ylim[2]-ylim[1])/10, paste("cvAUC = ",round(x$cvAUC,digits),sep=""), pos=4, cex=0.7, offset=0)
    }
  }
  )
  invisible()
}


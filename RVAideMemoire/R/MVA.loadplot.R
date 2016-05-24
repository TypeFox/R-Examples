# Loading plot
#  - 1 possible plot:
#     * PCA (dudi.pca[ade4],prcomp[stats],princomp[stats],pca[mixOmics],pca[labdsv],rda[vegan])
#     * sPCA (spca[mixOmics])
#	* IPCA (ipca[mixOmics])
#	* sIPCA (sipca[mixOmics])
#     * LDA (lda[MASS],discrimin[ade4])
#	* PLS-DA (plsda[mixOmics])
#	* sPLS-DA (splsda[mixOmics])
#	* Multilevel (s)PLS-DA (multilevel[mixOmics])
#	* CPPLS (mvr[pls])
#	* PLSR (mvr[pls],pls[mixOmics],plsR[plsRglm])
#	* sPLSR (spls[mixOmics])
#	* Multilevel (s)PLSR (multilevel[mixOmics]) # mixOmics > 5.0.4
#	* PLS-GLR (plsRglm[plsRglm])
#	* PCR (mvr[pls])
#	* CDA (discrimin.coa[ade4])
#	* Non Symmetric COA (dudi.nsc[ade4])
#  - Possibly separated for multiple factors:
#     * MCA (dudi.acm[ade4])
#     * Mix analysis (dudi.mix[ade4],dudi.hillsmith[ade4])
#  - Different sets of points in the same space:
#	* PCIA (procuste[ade4])
#  - Constrained/unconstrained spaces:
#     * 1 possible plot per space:
#	    ¤ RDA (pcaiv[ade4],pcaivortho[ade4],rda[vegan])
#  - Spaces from different data sets:
#	* 2 spaces (X and Y):
#         ¤ CIA (coinertia[ade4])
#	    ¤ rCCorA (rcc[mixOmics])
#	    ¤ 2B-PLS (pls[mixOmics])
#	    ¤ 2B-sPLS (spls[mixOmics])
#	    ¤ Multilevel 2B-(s)PLS (multilevel[mixOmics]) # mixOmics > 5.0.4
#	* >=2 spaces (including DA)
#	    ¤ rGCCA (rgcca[RGCCA],wrapper.rgcca[mixOmics])
#	    ¤ sGCCA (sgcca[RGCCA],wrapper.sgcca[mixOmics])


MVA.loadplot <- function(x,xax=1,yax=2,fac=NULL,set=c(12,1,2),space=1,map=TRUE,xlab=NULL,ylab=NULL,main=NULL,
  points=TRUE,ident=TRUE,links=TRUE,line=TRUE,labels=NULL,main.pos=c("bottomleft","topleft","bottomright","topright"),
  main.cex=1.3,legend=FALSE,legend.pos=c("topleft","topright","bottomleft","bottomright"),legend.title=NULL,
  legend.lab=NULL,pch=16,cex=1,col=1,lwd=1,lty=1,drawextaxes=TRUE,drawintaxes=TRUE,xlim=NULL,ylim=NULL) {
  lo <- MVA.load(x,xax,yax,set,space)
  loads <- lo$loads
  if (ncol(loads)<2) {map <- FALSE}
  if (is.null(labels)) {labels <- rownames(loads)}
  main.pos <- match.arg(main.pos)
  legend.pos <- match.arg(legend.pos)
  legend.col <- col
  legend.lwd <- lwd
  legend.lty <- lty
  if (map) {
    if (is.null(fac) & "set" %in% names(lo)) {fac <- lo$set}
    if (!is.null(fac)) {
	fac <- droplevels(factor(fac))
	if (is.null(legend.lab)) {legend.lab <- levels(fac)}
	if (length(legend.lab)!=nlevels(fac)) {stop("non-convenient 'legend.lab' argument")}
	if (length(pch)!=nlevels(fac)) {
	  if (length(pch)!=1) {stop("non-convenient 'pch' argument")}
	} else {
	  pch <- pch[as.numeric(fac)]
	}
	if (length(cex)!=nlevels(fac)) {
	  if (length(cex)!=1) {stop("non-convenient 'cex' argument")}
	} else {
	  cex <- cex[as.numeric(fac)]
	}
	if (length(col)!=nlevels(fac)) {
	  if (length(col)!=1) {stop("non-convenient 'col' argument")}
	} else {
	  col <- col[as.numeric(fac)]
	}
	if (length(lwd)!=nlevels(fac)) {
	  if (length(lwd)!=1) {stop("non-convenient 'lwd' argument")}
	} else {
	  lwd <- lwd[as.numeric(fac)]
	}
    }
  } else {
    if (ncol(loads)==2) {line <- TRUE}
    if (length(cex)>1) {cex <- cex[1]}
    if (ncol(loads)==2 & length(col)!=2) {
	if (length(col)==1) {
	  col <- rep(col,2)
	  legend.col <- col
	} else {stop("non-convenient 'col' argument")}
    }
    if (ncol(loads)==2 & length(lwd)!=2) {
	if (length(lwd)==1) {
	  lwd <- rep(lwd,2)
	  legend.lwd <- lwd
	} else {stop("non-convenient 'lwd' argument")}
    }
    if (ncol(loads)==2 & length(lty)!=2) {
	if (length(lty)==1) {
	  lty <- rep(lty,2)
	  legend.lty <- lty
	} else {stop("non-convenient 'lty' argument")}
    }
  }
  if (map) {
    MVA.loadplot.map(loads,xlab,ylab,main,points,ident,links,labels,main.pos,main.cex,legend,legend.pos,
	legend.title,legend.lab,legend.col,pch,cex,col,lwd,drawextaxes,drawintaxes,xlim,ylim)
  } else {
    MVA.loadplot.lin(loads,ylab,main,ident,line,labels,main.cex,legend,legend.pos,legend.title,legend.lab,
	legend.col,legend.lwd,legend.lty,cex,drawextaxes,ylim)
  }
}

MVA.loadplot.lin <- function(loads,ylab,main,ident,line,labels,main.cex,legend,
  legend.pos,legend.title,legend.lab,legend.col,legend.lwd,legend.lty,cex,drawextaxes,ylim) {
  oldmar <- par()$mar
  marinf <- ifelse(drawextaxes,5.1,3.5)
  if (drawextaxes) {
    par(mar=c(marinf,4.1,3.1,0.1))
  } else {
    par(mar=c(marinf,2.5,3.1,0.1))
  }
  if (is.null(ylim)) {ylim <- range(loads)}
  if (is.null(ylab)) {ylab <- "Loading"}
  plot(loads[,1],xlab="",ylab="",ylim=ylim,axes=FALSE,type="n")
  if(drawextaxes) {
    axis(2)
  }
  abline(h=0,col="grey",lty=3)
  lab.line <- ifelse(drawextaxes,2.3,0.8)
  mtext(ylab,side=2,line=lab.line,at=mean(ylim))
  if (ident) {
    axis(1,labels=labels,at=1:nrow(loads),las=2,cex.axis=cex)
  }
  if (ncol(loads)==1) {
    points(loads[,1],type=ifelse(line,"l","h"),col=legend.col,lwd=legend.lwd,lty=legend.lty)
  } else {
    for (i in 1:ncol(loads)) {lines(loads[,i],col=legend.col[i],lwd=legend.lwd[i],lty=legend.lty[i])}
  }
  if (!is.null(main)) {title(main,cex.main=main.cex)}
  if (legend) {
    if (is.null(legend.lab)) {legend.lab <- colnames(loads)}
    if (!is.null(legend.title) && nchar(legend.title)>0) {
	legend(legend.pos,legend.lab,col=legend.col,lwd=legend.lwd,lty=legend.lty,bg="white",title=legend.title)
    } else {
	legend(legend.pos,legend.lab,col=legend.col,lwd=legend.lwd,lty=legend.lty,bg="white")
    }
  }
  box()
  par(mar=oldmar)
}

MVA.loadplot.map <- function(loads,xlab,ylab,main,points,ident,links,labels,main.pos,main.cex,legend,
  legend.pos,legend.title,legend.lab,legend.col,pch,cex,col,lwd,drawextaxes,drawintaxes,xlim,ylim) {
  oldmar <- par()$mar
  marinf <- ifelse(drawextaxes,5.1,3.5)
  if (drawextaxes) {
    par(mar=c(marinf,4.1,2.1,0.1))
  } else {
    par(mar=c(marinf,2.5,2.1,0.1))
  }
  loadx <- loads[,1]
  loady <- loads[,2]
  if (is.null(xlab)) {xlab <- colnames(loads)[1]}
  if (is.null(ylab)) {ylab <- colnames(loads)[2]}
  if (is.null(xlim)) {
    xlim <- 1.1*range(loadx)
    if (all(xlim>0)) {xlim[1] <- -xlim[2]}
    if (all(xlim<0)) {xlim[2] <- -xlim[1]}
  }
  if (is.null(ylim)) {
    ylim <- 1.1*range(loady)
    if (all(ylim>0)) {ylim[1] <- -ylim[2]}
    if (all(ylim<0)) {ylim[2] <- -ylim[1]}
  }
  plot(loadx,loady,xlab="",ylab="",xlim=xlim,ylim=ylim,axes=FALSE,type="n")
  if(drawextaxes) {
    axis(1)
    axis(2)
  }
  if (drawintaxes) {abline(v=0,h=0,col="grey")}
  lab.line <- c(ifelse(drawextaxes,3,1),ifelse(drawextaxes,2.3,0.8))
  mtext(c(xlab,ylab),side=c(1,2),line=lab.line,at=c(mean(range(xlim)),mean(range(ylim))))
  if (points) {
    if (links) {
	segments(0,0,loadx,loady,col=col,lwd=lwd)
    } else {
	points(loadx,loady,pch=pch,cex=cex,col=col)
    }
    if (ident) {
	pos.lab <- function(x1,y1,x2,y2) {
	  res <- integer(length(x1))
	  for (i in 1:length(x1)) {
	    x2.i <- x2[i]-x1[i]
	    y2.i <- y2[i]-y1[i]
	    res[i] <- if (y2.i<0) {
		if (x2.i<0) {ifelse(atan(abs(y2.i/x2.i))*180/pi<40,2,1)} else {ifelse(atan(abs(y2.i/x2.i))*180/pi<40,4,1)}
	    } else {
		if (x2.i<0) {ifelse(atan(abs(y2.i/x2.i))*180/pi<40,2,3)} else {ifelse(atan(abs(y2.i/x2.i))*180/pi<40,4,3)}
	    }
	  }
	  return(res)
	}
	pos <- pos.lab(rep(0,length(loadx)),rep(0,length(loady)),loadx,loady)
	text(loadx,loady,labels,cex=cex,col=col,pos=pos,offset=0.2)
    }
  } else {
    text(loadx,loady,labels,cex=cex,col=col)
  }
  if (!is.null(main)) {
    xmain <- if (main.pos %in% c("bottomleft","topleft")) {xlim[1]-0.02*diff(xlim)} else {xlim[2]+0.02*diff(xlim)}
    ymain <- if (main.pos %in% c("bottomleft","bottomright")) {ylim[1]} else {ylim[2]}
    adjmain <- if (main.pos %in% c("bottomleft","topleft")) {c(0,NA)} else {c(1,NA)}
    text(xmain,ymain,main,adj=adjmain,cex=main.cex)
  }
  if (legend) {
    if (is.null(legend.lab)) {legend.lab <- "1"}
    if (!is.null(legend.title) && nchar(legend.title)>0) {
	legend(legend.pos,legend.lab,fill=legend.col,bg="white",title=legend.title)
    } else {
	legend(legend.pos,legend.lab,fill=legend.col,bg="white")
    }
  }
  box()
  par(mar=oldmar)
}




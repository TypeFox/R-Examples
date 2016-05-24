# Score plot
#  - 1 possible plot:
#     * PCA (dudi.pca[ade4],prcomp[stats],princomp[stats]°,pca[mixOmics],pca[labdsv],rda[vegan])
#		° if scores=TRUE
#     * sPCA (spca[mixOmics])
#	* IPCA (ipca[mixOmics])
#	* sIPCA (sipca[mixOmics])
#     * PCoA (dudi.pco[ade4],pcoa[ape],pco[labdsv],cmdscale[stats]°,wcmdscale[vegan]°,capscale[vegan])
#		° if computed with at least one non-default argument
#	* nMDS (isoMDS[MASS],monoMDS[vegan],metaMDS[vegan],nmds[labdsv])
#     * LDA (lda[MASS],discrimin[ade4])
#	* PLS-DA (plsda[mixOmics])
#	* sPLS-DA (splsda[mixOmics])
#	* Multilevel (s)PLS-DA (multilevel[mixOmics]) # mixOmics >= 5.0.4
#	* CPPLS (mvr[pls])
#	* PLSR (mvr[pls],pls[mixOmics],plsR[plsRglm])
#	* sPLSR (spls[mixOmics])
#	* Multilevel (s)PLSR (multilevel[mixOmics]) # mixOmics >= 5.0.4
#	* PLS-GLR (plsRglm[plsRglm])
#	* PCR (mvr[pls])
#	* CDA (discrimin[ade4],discrimin.coa[ade4])
#	* Non Symmetric COA (dudi.nsc[ade4])
#  - Possibly separated for multiple factors:
#     * MCA (dudi.acm[ade4])
#     * Mix analysis (dudi.mix[ade4],dudi.hillsmith[ade4])
#  - Different sets of points in the same space:
#	* COA (dudi.coa[ade4],cca[vegan])
#	* Decentred COA (dudi.dec[ade4])
#	* PCIA (procuste[ade4])
#     * DPCoA (dpcoa[ade4])
#  - Constrained/unconstrained spaces:
#     * 1 possible plot per space:
#	    ¤ RDA (pcaiv[ade4],pcaivortho[ade4],rda[vegan])
#	    ¤ db-RDA (capscale[vegan])
#     * Rows and/or columns per space:
#         ¤ CCA (cca[vegan],cca[ade4])
#  - Spaces from different data sets:
#	* 2 spaces (X and Y):
#         ¤ CCorA (CCorA[vegan])
#	* 3 spaces (X, Y and "common"):
#	    ¤ CCorA (rcc[mixOmics])
#	    ¤ rCCorA (rcc[mixOmics])
#         ¤ CIA (coinertia[ade4])
#	    ¤ 2B-PLS (pls[mixOmics])
#	    ¤ 2B-sPLS (spls[mixOmics])
#	    ¤ Multilevel 2B-(s)PLS (multilevel[mixOmics]) # mixOmics >= 5.0.4
#	* >=2 spaces (including DA)
#	    ¤ rGCCA (rgcca[RGCCA],wrapper.rgcca[mixOmics])
#	    ¤ sGCCA (sgcca[RGCCA],wrapper.sgcca[mixOmics])


MVA.scoreplot <- function(x,xax=1,yax=2,scaling=2,set=c(12,1,2),space=1,byfac=TRUE,fac=NULL,barycenters=TRUE,
  stars=TRUE,contours=FALSE,dhist=TRUE,weights=1,xlab=NULL,ylab=NULL,main=NULL,pch=16,cex=1,col=1,points=TRUE,
  labels=NULL,main.pos=c("bottomleft","topleft","bottomright","topright"),main.cex=1.3,fac.lab=NULL,fac.cex=1,
  legend=FALSE,legend.pos=c("topleft","topright","bottomleft","bottomright"),legend.title=NULL,legend.lab=NULL,
  legend.cex=1,drawextaxes=TRUE,drawintaxes=TRUE,xlim=NULL,ylim=NULL,keepmar=FALSE) {
  sco <- MVA.scores(x,xax,yax,scaling,set,space)
  coord <- sco$coord
  if (length(weights)!=nrow(coord)) {
    if (length(weights)==1) {
	weights=rep(weights,nrow(coord))
    } else {stop("non-convenient 'weights' argument")}
  }
  if ((inherits(x,c("COA.ade4","DCOA.ade4","CCA.ade4","COA.vegan","CCA.vegan","PCIA.ade4","DPCoA.ade4")) && set==12) |
    (inherits(x,"CIA.ade4") && space==3 && set==12)) {
    fac <- sco$set
    barycenters <- contours <- stars <- FALSE
    if (length(cex)==1) {cex <- c(1,0.7)}
    if (length(pch)==1) {pch <- c(16,17)}
  }
  if (!points & is.null(labels)) {labels <- rownames(coord)}
  main.pos <- match.arg(main.pos)
  legend.pos <- match.arg(legend.pos)
  legend.col <- col
  legend.pch <- pch
  legend.ptcex <- cex
  if (!is.null(fac)) {
    fac <- droplevels(factor(fac))
    if (is.null(fac.lab)) {fac.lab <- levels(fac)}
    if (!length(fac.cex) %in% c(1,nlevels(fac))) {stop("non-convenient 'fac.cex' argument")}
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
  }
  if (inherits(x,c("MCA.ade4","Mix.ade4")) && byfac) {
    MVA.scoreplot.byfac(coord,x,barycenters,stars,contours,weights,xlab,ylab,main,pch,cex,col,
	points,labels,main.pos,main.cex,fac.lab,fac.cex,legend,legend.pos,legend.title,legend.lab,legend.col,
	legend.pch,legend.ptcex,legend.cex,drawintaxes,xlim,ylim,keepmar)
  } else if (inherits(x,c("COA.ade4","DCOA.ade4","CCA.ade4","COA.vegan","CCA.vegan","PCIA.ade4","DPCoA.ade4")) && ncol(coord)==1) {
    if (dhist) {
	MVA.scoreplot.1comp.dhist(coord,fac,xlab,ylab,legend.col,legend,legend.pos,legend.title,legend.lab,
	  drawextaxes,drawintaxes,main,xlim)
    } else {
	MVA.scoreplot.1comp.dotchart(coord,fac,barycenters,weights,xlab,main,pch,cex,col,labels,legend.col,
	  drawintaxes,xlim)
    }
  } else {
    if (ncol(coord)==1) {
	if (dhist) {
	  MVA.scoreplot.1comp.dhist(coord,fac,xlab,ylab,legend.col,legend,legend.pos,legend.title,legend.lab,
	    drawextaxes,drawintaxes,main,xlim)
	} else {
	  MVA.scoreplot.1comp.dotchart(coord,fac,barycenters,weights,xlab,main,pch,cex,col,labels,legend.col,
	    drawintaxes,xlim)
	}
    } else {
	MVA.scoreplot.2comp(coord,byfac=FALSE,fac,barycenters,stars,contours,weights,xlab,ylab,main,pch,cex,col,
	  points,labels,main.pos,main.cex,fac.lab,fac.cex,legend,legend.pos,legend.title,legend.lab,legend.col,
	  legend.pch,legend.ptcex,legend.cex,drawextaxes,drawintaxes,xlim,ylim,keepmar)
    }
  }
}

MVA.scoreplot.byfac <- function(coord,x,barycenters,stars,contours,weights,xlab,ylab,main,pch,cex,col,points,
  labels,main.pos,main.cex,fac.lab,fac.cex,legend,legend.pos,legend.title,legend.lab,legend.col,legend.pch,
  legend.ptcex,legend.cex,drawintaxes,xlim,ylim,keepmar) {
  oritab <- eval.parent(as.list(x$call)[[2]])
  if (inherits(MVA.ident(x),"Mix.ade4")) {
    varnames <- colnames(oritab)
    if (length(which(x$index %in% c("f","o")))==0) {stop("no factor in the analysis")}
    oritab <- oritab[,which(x$index %in% c("f","o"))]
    if (!is.data.frame(oritab)) {
	oritab <- as.data.frame(oritab)
	colnames(oritab) <- varnames[which(x$index %in% c("f","o"))]
    }
  }
  nf <- ncol(oritab)
  if (length(legend.col)==1) {legend.col <- rep(legend.col[1],max(apply(oritab,2,function(y) nlevels(factor(y)))))}
  oldmfrow <- par()$mfrow
  par(mfrow=n2mfrow(nf))
  for (i in 1:nf) {
    fac <- factor(oritab[,i])
    if (ncol(coord)==1) {
	MVA.scoreplot.1comp.dhist(coord,fac,xlab,ylab,col=legend.col[1:nlevels(fac)],legend,legend.pos,legend.title=NULL,
	  legend.lab=levels(fac),drawextaxes=FALSE,drawintaxes,main=ifelse(is.null(main),colnames(oritab)[i],main[i]),xlim)
    } else {
	MVA.scoreplot.2comp(coord,byfac=TRUE,fac,barycenters,stars,contours,weights,xlab,ylab,main=ifelse(is.null(main),
	  colnames(oritab)[i],main[i]),pch,cex,col=legend.col[1:nlevels(fac)][as.numeric(fac)],points,labels,
	  main.pos,main.cex,fac.lab=levels(fac),fac.cex,legend,legend.pos,legend.title=NULL,legend.lab=levels(fac),
	  legend.col=legend.col[1:nlevels(fac)],legend.pch,legend.ptcex,legend.cex,drawextaxes=FALSE,drawintaxes,
	  xlim,ylim,keepmar)
    }
  }
  par(mfrow=oldmfrow)
}

MVA.scoreplot.1comp.dhist <- function(coord,fac,xlab,ylab,col,legend,legend.pos,legend.title,legend.lab,drawextaxes,
  drawintaxes,main,xlim) {
  if (is.null(fac)) {
    fac <- gl(1,nrow(coord))
  } else {
    if (length(col)==1) {col <- rep(col,nlevels(fac))}
  }
  if (is.null(xlab)) {xlab <- colnames(coord)}
  if (is.null(xlim)) {xlim <- range(coord)}
  dhist(as.vector(t(as.matrix(coord))),fac=fac,col=col,legend=legend,pos.legend=legend.pos,title.legend=legend.title,
    lab.legend=legend.lab,xlab=xlab,ylab=ylab,drawextaxes=drawextaxes,drawintaxes=drawintaxes,main=main,xlim=xlim)
}

MVA.scoreplot.1comp.dotchart <- function(coord,fac,barycenters,weights,xlab,main,pch,cex,col,labels,legend.col,
  drawintaxes,xlim) {
  coordx <- coord[,1]
  names(coordx) <- rownames(coord)
  if (is.null(xlab)) {xlab <- colnames(coord)}
  if (barycenters) {
    bary <- numeric(nlevels(fac))
    for (i in 1:nlevels(fac)) {
	bary[i] <- wmean(coordx[as.numeric(fac)==i],weights[as.numeric(fac)==i])
    }
  } else {bary <- NULL}
  if (length(cex)>1) {cex <- cex[1]}
  if (is.null(xlim)) {xlim <- range(coordx)}
  dotchart(coordx,labels=labels,groups=fac,gdata=bary,cex=cex,pch=pch,gpch=15,color=col,gcolor=legend.col,
    main=main,xlab=xlab,xlim=xlim)
  if (drawintaxes) {abline(v=0,col="grey")}
}

MVA.scoreplot.2comp <- function(coord,byfac,fac,barycenters,stars,contours,weights,xlab,ylab,main,pch,cex,col,
  points,labels,main.pos,main.cex,fac.lab,fac.cex,legend,legend.pos,legend.title,legend.lab,legend.col,legend.pch,
  legend.ptcex,legend.cex,drawextaxes,drawintaxes,xlim,ylim,keepmar) {
  oldmar <- par()$mar
  marsup <- ifelse(byfac,0,2.1)
  marinf <- if (byfac) {2.5} else {ifelse(drawextaxes,5.1,3.5)}
  if (drawextaxes) {
    par(mar=c(marinf,4.1,marsup,0.1))
  } else {
    par(mar=c(marinf,2.5,marsup,0.1))
  }
  coordx <- coord[,1]
  coordy <- coord[,2]
  if (is.null(xlab)) {xlab <- colnames(coord)[1]}
  if (is.null(ylab)) {ylab <- colnames(coord)[2]}
  if (is.null(xlim)) {xlim <- range(coordx)}
  if (is.null(ylim)) {ylim <- range(coordy)}
  plot(coordx,coordy,xlab="",ylab="",xlim=xlim,ylim=ylim,axes=FALSE,type="n")
  if(drawextaxes) {
    axis(1)
    axis(2)
  }
  if (drawintaxes) {abline(v=0,h=0,col="grey")}
  lab.line <- c(ifelse(drawextaxes,3,1),ifelse(drawextaxes,2.3,0.8))
  mtext(c(xlab,ylab),side=c(1,2),line=lab.line,at=c(mean(range(coordx)),mean(range(coordy))))
  if (!is.null(fac)) {
    fac.col <- legend.col
    if (length(fac.col)==1) {fac.col <- rep(fac.col,nlevels(fac))}
    if (length(fac.cex)==1) {fac.cex <- rep(fac.cex,nlevels(fac))}
    bar.x <- bar.y <- numeric(nlevels(fac))
    for (i in 1:nlevels(fac)) {
	bar.x[i] <- wmean(coordx[as.numeric(fac)==i],weights[as.numeric(fac)==i])
	bar.y[i] <- wmean(coordy[as.numeric(fac)==i],weights[as.numeric(fac)==i])
    }
    if (contours) {
	for (i in 1:nlevels(fac)) {
	  coordx.temp <- coordx[as.numeric(fac)==i]
	  coordy.temp <- coordy[as.numeric(fac)==i]
	  long <- length(coordx.temp)
	  longinit <- long
	  cref <- 1
	  repeat {
	    if (long<3) {break}
	    if (cref==0) {break}
	    num <- chull(coordx.temp,coordy.temp)
	    x2 <- coordx.temp[num]
	    y2 <- coordy.temp[num]
	    taux <- long/longinit
	    if (taux<=cref & cref==1) {
		polygon(x2,y2,border=fac.col[i])
		cref <- 0
	    }
	    coordx.temp <- coordx.temp[-num]
	    coordy.temp <- coordy.temp[-num]
	    long <- length(coordx.temp)
	  }
	}
   }
    if (stars) {
	for (i in 1:nlevels(fac)) {
	  coordx.temp <- coordx[as.numeric(fac)==i]
	  coordy.temp <- coordy[as.numeric(fac)==i]
	  segments(bar.x[i],bar.y[i],coordx.temp,coordy.temp,col=fac.col[i])
	}
    }
  }
  if (points) {
    points(coordx,coordy,pch=pch,cex=cex,col=col)
  } else {
    text(coordx,coordy,labels,cex=cex,col=col)
  }
  if (!is.null(fac)) {
    if (barycenters) {
	lab <- paste0(" ",fac.lab," ")
	for (i in 1:nlevels(fac)) {
	  xh <- strwidth(lab[i],cex=fac.cex[i])
	  yh <- strheight(lab[i],cex=fac.cex[i])*5/3
	  rect(bar.x[i]-xh/2,bar.y[i]-yh/2,bar.x[i]+xh/2,bar.y[i]+yh/2,col="white",border=fac.col[i])
	  text(bar.x[i],bar.y[i],lab[i],col=fac.col[i],cex=fac.cex[i])
	}
    }
  }
  if (!is.null(main)) {
    xmain <- if (byfac) {
	if (main.pos %in% c("bottomleft","topleft")) {xlim[1]} else {xlim[2]}
    } else  {
	if (main.pos %in% c("bottomleft","topleft")) {xlim[1]-0.02*diff(xlim)} else {xlim[2]+0.02*diff(xlim)}
    }
    ymain <- if (byfac) {
	if (main.pos %in% c("bottomleft","bottomright")) {ylim[1]+0.05*diff(ylim)} else {ylim[2]-0.05*diff(ylim)}
    } else {
	if (main.pos %in% c("bottomleft","bottomright")) {ylim[1]} else {ylim[2]}
    }
    adjmain <- if (main.pos %in% c("bottomleft","topleft")) {c(0,NA)} else {c(1,NA)}
    text(xmain,ymain,main,adj=adjmain,cex=main.cex)
  }
  if (legend) {
    if (points) {
	if (!is.null(legend.title) && nchar(legend.title)>0) {
	  legend(legend.pos,legend.lab,col=legend.col,pch=legend.pch,pt.cex=legend.ptcex,
	    lty=ifelse(legend.ptcex==0,1,0),cex=legend.cex,bg="white",title=legend.title,title.col=1)
	} else {
	  legend(legend.pos,legend.lab,col=legend.col,pch=legend.pch,pt.cex=legend.ptcex,
	    lty=ifelse(legend.ptcex==0,1,0),cex=legend.cex,bg="white")
	}
    } else {
	if (!is.null(legend.title) && nchar(legend.title)>0) {
	  legend(legend.pos,legend.lab,text.col=legend.col,cex=legend.cex,bg="white",title=legend.title,title.col=1)
	} else {
	  legend(legend.pos,legend.lab,text.col=legend.col,cex=legend.cex,bg="white")
	}
    }
  }
  box()
  if (!keepmar) {par(mar=oldmar)}
  invisible(list(xlim=xlim,ylim=ylim))
}

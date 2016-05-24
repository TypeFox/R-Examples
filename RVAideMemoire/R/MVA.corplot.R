# Correlation circle
#  - 1 possible plot:
#     * PCA (dudi.pca[ade4],pca[mixOmics],rda[vegan])
#     * sPCA (spca[mixOmics])
#     * IPCA (ipca[mixOmics])
#     * sIPCA (sipca[mixOmics])
#	* LDA (lda[MASS],discrimin[ade4])
#	* PLS-DA (plsda[mixOmics])
#	* sPLS-DA (splsda[mixOmics])
#	* Multilevel (s)PLS-DA (multilevel[mixOmics]) # mixOmics >= 5.0.4
#	* CPPLS (mvr[pls])
#	* PLSR (mvr[pls],pls[mixOmics],plsR[plsRglm])
#	* sPLSR (spls[mixOmics])
#	* Multilevel (s)PLSR (multilevel[mixOmics]) # mixOmics > 5.0.4
#	* PLS-GLR (plsRglm[plsRglm])
#	* PCR (mvr[pls])
#	* CDA (discrimin.coa[ade4])
#	* Non Symmetric COA (dudi.nsc[ade4]) # Pas corrélation
#	* CCA (cca[vegan],cca[ade4])
#  - Possibly separated for multiple factors:
#     * Mix analysis (dudi.mix[ade4],dudi.hillsmith[ade4])
#  - Constrained/unconstrained spaces:
#     * 1 possible plot per space:
#	    ¤ RDA (rda[vegan],pcaiv[ade4],pcaivortho[ade4])
#  - Spaces from different data sets:
#	* 2 spaces (X and Y):
#         ¤ CCorA (CCorA[vegan],rcc[mixOmics])
#	* 3 spaces (X, Y and "common"):
#         ¤ CIA (coinertia[ade4])
#	    ¤ rCCorA (rcc[mixOmics])
#	    ¤ 2B-PLS (pls[mixOmics])
#	    ¤ 2B-sPLS (spls[mixOmics])
#	    ¤ Multilevel 2B-(s)PLS (multilevel[mixOmics]) # mixOmics > 5.0.4
#	* >=2 spaces (including DA)
#	    ¤ rGCCA (wrapper.rgcca[mixOmics])
#	    ¤ sGCCA (wrapper.sgcca[mixOmics])


MVA.corplot <- function(x,xax=1,yax=2,thresh=0,fac=NULL,set=c(12,1,2),space=1,xlab=NULL,ylab=NULL,main=NULL,
  circle=TRUE,intcircle=0.5,points=TRUE,ident=TRUE,arrows=TRUE,labels=NULL,main.pos=c("bottomleft","topleft",
  "bottomright","topright"),main.cex=1.3,legend=FALSE,legend.pos=c("topleft","topright","bottomleft","bottomright"),
  legend.title=NULL,legend.lab=NULL,pch=16,cex=1,col=1,lwd=1,drawintaxes=TRUE,add=FALSE,add.const=1,keepmar=FALSE) {
  co <- MVA.cor(x,xax,yax,set,space)
  corr <- co$corr
  if (is.numeric(thresh) && thresh!=0) {
    to.keep <- apply(corr,1,function(x) any(abs(x)>thresh))
    if (sum(to.keep)==0) {stop("too stringent threshold, no variable retained")}
    corr <- corr[to.keep,]
    if ("set" %in% names(co)) {co$set <- co$set[to.keep]}
  }
  if (is.null(labels)) {labels <- rownames(corr)}
  main.pos <- match.arg(main.pos)
  legend.pos <- match.arg(legend.pos)
  legend.col <- col
  if (is.null(fac) & "set" %in% names(co)) {fac <- co$set}
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
  if (ncol(corr)==1) {
    MVA.corplot.1comp(corr,fac,xlab,arrows,main,pch,cex,col,lwd,labels,legend.col,drawintaxes)
  } else {
    MVA.corplot.2comp(corr,xlab,ylab,main,circle,intcircle,points,ident,arrows,labels,main.pos,main.cex,
	legend,legend.pos,legend.title,legend.lab,legend.col,pch,cex,col,lwd,drawintaxes,add,add.const,
	keepmar)
  }
}

MVA.corplot.1comp <- function(corr,fac,xlab,arrows,main,pch,cex,col,lwd,labels,legend.col,drawintaxes) {
  corrx <- corr[,1]
  names(corrx) <- rownames(corr)
  if (is.null(xlab)) {xlab <- colnames(corr)}
  if (length(cex)>1) {cex <- cex[1]}
  lim <- max(abs(corrx))
  if (lim<1) {lim <- 1}
  dotchart(corrx,labels=labels,groups=fac,cex=cex,pch=NA,gpch=15,color=col,gcolor=legend.col,
    main=main,xlab=xlab,xlim=c(-lim,lim))
  if (drawintaxes) {abline(v=0,col="grey")}
  if (arrows) {
    if (!is.null(fac)) {
	y <- length(corrx)+2*(nlevels(fac)-1)
	for (i in 1:nlevels(fac)) {
	  col.temp <- if (length(unique(col))==nlevels(fac)) {unique(col[as.numeric(fac)==i])} else {col}
	  lwd.temp <- if (length(unique(lwd))==nlevels(fac)) {unique(lwd[as.numeric(fac)==i])} else {lwd}
	  arrows(0,y:(y+1-length(corrx[as.numeric(fac)==i])),rev(corrx[as.numeric(fac)==i]),
	    y:(y+1-length(corrx[as.numeric(fac)==i])),length=0.08,col=col.temp,lwd=lwd.temp)
	  y <- y-length(corrx[as.numeric(fac)==i])-2
	}
    } else {
	arrows(0,1:length(corrx),corrx,1:length(corrx),length=0.08,col=col,lwd=lwd)
    }
  } else {
    if (!is.null(fac)) {
	y <- length(corrx)+2*(nlevels(fac)-1)
	for (i in 1:nlevels(fac)) {
	  pch.temp <- if (length(unique(pch))==nlevels(fac)) {unique(pch[as.numeric(fac)==i])} else {pch}
	  col.temp <- if (length(unique(col))==nlevels(fac)) {unique(col[as.numeric(fac)==i])} else {col}
	  points(rev(corrx[as.numeric(fac)==i]),y:(y+1-length(corrx[as.numeric(fac)==i])),pch=pch.temp,cex=cex,
	    col=col.temp)
	  y <- y-length(corrx[as.numeric(fac)==i])-2
	}
    } else {
	points(corrx,1:length(corrx),pch=pch,cex=cex,col=col)
    }
  }
}

MVA.corplot.2comp <- function(corr,xlab,ylab,main,circle,intcircle,points,ident,arrows,labels,main.pos,
  main.cex,legend,legend.pos,legend.title,legend.lab,legend.col,pch,cex,col,lwd,drawintaxes,add,add.const,
  keepmar) {
  if (add) {corr <- corr*add.const}
  corrx <- corr[,1]
  corry <- corr[,2]
  oldmar <- par()$mar
  lim <- max(abs(corr))
  if (lim<1) {lim <- 1}
  if (!add) {
    par(mar=c(3,3,1,1))
    plot(corrx,corry,xlim=c(-lim,lim),ylim=c(-lim,lim),type="n",axes=FALSE,xlab="",ylab="",xaxs="i",
	yaxs="i")
  } else {
    main <- NULL
  }
  if (is.null(xlab)) {xlab <- colnames(corr)[1]}
  if (is.null(ylab)) {ylab <- colnames(corr)[2]}
  mtext(c(xlab,ylab),side=c(1,2),line=1.2,at=0)
  if (drawintaxes) {abline(v=0,h=0,col="grey")}
  if (points) {
    if (arrows) {
	arrows(0,0,corrx,corry,length=0.08,col=col,lwd=lwd)
    } else {
	points(corrx,corry,pch=pch,cex=cex,col=col)
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
	pos <- pos.lab(rep(0,length(corrx)),rep(0,length(corry)),corrx,corry)
	text(corrx,corry,labels,cex=cex,col=col,pos=pos,offset=0.2,xpd=TRUE)
    }
  } else {
    text(corrx,corry,labels,cex=cex,col=col,xpd=TRUE)
  }
  if (circle) {symbols(0,0,circles=1,inches=FALSE,add=TRUE)}
  if (!is.null(intcircle) && !is.na(intcircle) && !intcircle %in% c(0,1)) {
    for (i in 1:length(intcircle)) {
	symbols(0,0,circles=intcircle[i],inches=FALSE,lty=3,add=TRUE)
    }
  }
  if (!is.null(main)) {
    xmain <- if (main.pos %in% c("bottomleft","topleft")) {-lim} else {lim}
    ymain <- if (main.pos %in% c("bottomleft","bottomright")) {-lim} else {lim-0.01*2*lim}
    adjmain <- if (main.pos %in% c("bottomleft","topleft")) {c(0,NA)} else {c(1,NA)}
    text(xmain,ymain,main,adj=adjmain,cex=main.cex,xpd=TRUE)
  }
  if (legend) {
    if (is.null(legend.lab)) {legend.lab <- "1"}
    if (!add) {
	par(new=TRUE)
	par(mar=c(0,0,0,0))
	plot(1,type="n",xlab="",ylab="",axes=FALSE)
    }
    if (!is.null(legend.title) && nchar(legend.title)>0) {
	legend(legend.pos,legend.lab,fill=legend.col,bg="white",title=legend.title)
    } else {
	legend(legend.pos,legend.lab,fill=legend.col,bg="white")
    }
  }
  if (!keepmar) {par(mar=oldmar)}
}


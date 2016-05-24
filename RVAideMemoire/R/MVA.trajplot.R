MVA.trajplot <- function(x,xax=1,yax=2,trajects,trajlab=NULL,scaling=2,set=c(12,1,2),space=1,xlab=NULL,ylab=NULL,
  main=NULL,pch=16,cex=1,trajlab.cex=1,col=1,lwd=1,lty=1,points=TRUE,allpoints=TRUE,arrows=TRUE,labels=NULL,
  main.pos=c("bottomleft","topleft","bottomright","topright"),main.cex=1.3,legend=FALSE,legend.pos=c("topleft",
  "topright","bottomleft","bottomright"),legend.title=NULL,legend.lab=NULL,legend.cex=1,drawextaxes=TRUE,
  drawintaxes=TRUE,xlim=NULL,ylim=NULL) {
  if (!is.list(trajects)) {trajects <- list(trajects)}
  if (ncol(as.data.frame(trajects[[1]]))!=1) {stop("wrong 'trajects', must be a vector or a list of vectors")}
  coord <- MVA.scores(x,xax,yax,scaling,set,space)$coord
  if (ncol(coord)==1) {stop("choose a second axis")}
  ntraj <- length(trajects)
  if (!is.null(trajlab) & length(trajlab)!=ntraj) {stop("non-convenient 'trajlab' argument")}
  which.in <- unique(unlist(trajects))
  rest <- !length(which.in)==nrow(coord)
  if (length(col)==1) {
    if (rest) {col <- rep(col,ntraj+1)} else {col <- rep(col,ntraj)}
  }
  if ((rest & length(col)!=ntraj+1) | (!rest & length(col)!=ntraj)) {stop("non-convenient 'col' argument")}
  if (length(lwd)==1) {
    if (rest) {lwd <- rep(lwd,ntraj+1)} else {lwd <- rep(lwd,ntraj)}
  }
  if ((rest & length(lwd)!=ntraj+1) | (!rest & length(lwd)!=ntraj)) {stop("non-convenient 'lwd' argument")}
  if (length(lty)==1) {
    if (rest) {lty <- rep(lty,ntraj+1)} else {lty <- rep(lty,ntraj)}
  }
  if ((rest & length(lty)!=ntraj+1) | (!rest & length(lty)!=ntraj)) {stop("non-convenient 'lty' argument")}
  if (length(trajlab.cex)==1) {trajlab.cex <- rep(trajlab.cex,ntraj)}
  if (length(trajlab.cex)!=ntraj) {stop("non-convenient 'trajlab.cex' argument")}
  if (points) {
    if (length(pch)==1) {
	if (rest) {pch <- rep(pch,ntraj+1)} else {pch <- rep(pch,ntraj)}
    }
    if ((rest & length(pch)!=ntraj+1) | (!rest & length(pch)!=ntraj)) {stop("non-convenient 'pch' argument")}
  }
  if (!points & is.null(labels)) {labels <- rownames(coord)}
  main.pos <- match.arg(main.pos)
  legend.pos <- match.arg(legend.pos)
  oldmar <- par()$mar
  marinf <- ifelse(drawextaxes,5.1,3.5)
  if (drawextaxes) {
    par(mar=c(marinf,4.1,2.1,0.1))
  } else {
    par(mar=c(marinf,2.5,2.1,0.1))
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
  if (rest & allpoints) {
    if (points) {
	points(coordx[-which.in],coordy[-which.in],pch=pch[ntraj+1],col=col[ntraj+1])
    } else {
	text(coordx[-which.in],coordy[-which.in],labels[-which.in],col=col[ntraj+1],cex=cex)
    }
  }
  for (i in 1:ntraj) {
    traj.i <- trajects[[i]]
    n <- length(traj.i)
    if (arrows) {
	cx <- coordx[traj.i]
	cy <- coordy[traj.i]
	arrows(cx[-n],cy[-n],cx[-n]+diff(cx)/2,cy[-n]+diff(cy)/2,col=col[i],lwd=lwd[i],lty=lty[i],
	  length=0.12,angle=20)
	segments(cx[-n]+diff(cx)/2,cy[-n]+diff(cy)/2,coordx[traj.i[-1]],coordy[traj.i[-1]],col=col[i],
	  lwd=lwd[i],lty=lty[i])
    } else {
	segments(coordx[traj.i[-n]],coordy[traj.i[-n]],coordx[traj.i[-1]],coordy[traj.i[-1]],col=col[i],
	  lwd=lwd[i],lty=lty[i])
    }
    if (points) {
	points(coordx[traj.i],coordy[traj.i],pch=pch[i],col=col[i])
    } else {
	text(coordx[traj.i],coordy[traj.i],labels[traj.i],col=col[i],cex=cex)
    }
    if (!is.null(trajlab)) {
	tlx <- cx[-n]+diff(cx)/2
	tly <- cy[-n]+diff(cy)/2
	wh <- if (length(tlx)%%2==0) {length(tlx)/2} else {(length(tlx)+1)/2}
	ctlx <- tlx[wh]
	ctly <- tly[wh]
	lab <- paste0(" ",trajlab[i]," ")
	xh <- strwidth(lab,cex=trajlab.cex[i])
	yh <- strheight(lab,cex=trajlab.cex[i])*5/3
	rect(ctlx-xh/2,ctly-yh/2,ctlx+xh/2,ctly+yh/2,col="white",border=col[i])
	text(ctlx,ctly,lab,col=col[i],cex=trajlab.cex[i])
    }
  }
  if (!is.null(main)) {
    xmain <- if (main.pos %in% c("bottomleft","topleft")) {xlim[1]-0.02*diff(xlim)} else {xlim[2]+0.02*diff(xlim)}
    ymain <- if (main.pos %in% c("bottomleft","bottomright")) {ylim[1]} else {ylim[2]}
    adjmain <- if (main.pos %in% c("bottomleft","topleft")) {c(0,NA)} else {c(1,NA)}
    text(xmain,ymain,main,adj=adjmain,cex=main.cex)
  }
  if (legend) {
    if (is.null(legend.lab)) {
	if (!is.null(trajlab)) {legend.lab <- trajlab} else {legend.lab <- as.character(1:ntraj)}
    }
    if (points) {
	if (!is.null(legend.title) && nchar(legend.title)>0) {
	  legend(legend.pos,legend.lab,col=col,pch=pch,lty=lty,cex=legend.cex,bg="white",title=legend.title)
	} else {
	  legend(legend.pos,legend.lab,col=col,pch=pch,lty=lty,cex=legend.cex,bg="white")
	}
    } else {
	if (!is.null(legend.title) && nchar(legend.title)>0) {
	  legend(legend.pos,legend.lab,col=col,lty=lty,cex=legend.cex,bg="white",title=legend.title)
	} else {
	  legend(legend.pos,legend.lab,col=col,lty=lty,cex=legend.cex,bg="white")
	}
    }
  }
  box()
  par(mar=oldmar)
}

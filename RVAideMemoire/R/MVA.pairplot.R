MVA.pairplot <- function(x,xax=1,yax=2,pairs=NULL,scaling=2,space=1,fac=NULL,xlab=NULL,ylab=NULL,main=NULL,
  ident=TRUE,labels=NULL,cex=0.7,col=1,lwd=1,main.pos=c("bottomleft","topleft","bottomright","topright"),
  main.cex=1.3,legend=FALSE,legend.pos=c("topleft","topright","bottomleft","bottomright"),legend.title=NULL,
  legend.lab=NULL,drawextaxes=TRUE,drawintaxes=TRUE,xlim=NULL,ylim=NULL) {
  sco <- MVA.scores(x,xax,yax,scaling,set=12,space)
  coord <- sco$coord
  if (ncol(coord)==1) {stop("choose a second axis")}
  if (is.null(pairs)) {
    if (!"set" %in% names(sco)) {
	stop("unknown relationships between points, use 'pairs' argument")
    } else {
	pairs <- sco$set
    }
  }
  if (!is.factor(pairs)) {pairs <- factor(pairs)}
  if (nlevels(pairs)!=2) {stop("there has to be two sets of points")}
  if (diff(table(pairs))!=0 | length(pairs)!=nrow(coord)) {stop("unclear relationships between points")}
  legend.pos <- match.arg(legend.pos)
  legend.col <- col
  legend.lwd <- lwd
  if (!is.null(fac)) {
    fac <- droplevels(factor(fac))
    if (is.null(legend.lab)) {legend.lab <- levels(fac)}
    if (length(legend.lab)!=nlevels(fac)) {stop("non-convenient 'legend.lab' argument")}
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
  oldmar <- par()$mar
  marinf <- ifelse(drawextaxes,5.1,3.5)
  if (drawextaxes) {
    par(mar=c(marinf,4.1,2.1,0.1))
  } else {
    par(mar=c(marinf,2.5,2.1,0.1))
  }
  coordx1 <- coord[,1][as.numeric(pairs)==1]
  coordy1 <- coord[,2][as.numeric(pairs)==1]
  coordx2 <- coord[,1][as.numeric(pairs)==2]
  coordy2 <- coord[,2][as.numeric(pairs)==2]
  if (is.null(xlab)) {xlab <- colnames(coord)[1]}
  if (is.null(ylab)) {ylab <- colnames(coord)[2]}
  if (is.null(xlim)) {xlim <- range(coord[,1])}
  if (is.null(ylim)) {ylim <- range(coord[,2])}
  plot(coord[,1],coord[,2],xlab="",ylab="",xlim=xlim,ylim=ylim,axes=FALSE,type="n")
  if(drawextaxes) {
    axis(1)
    axis(2)
  }
  if (drawintaxes) {abline(v=0,h=0,col="grey")}
  lab.line <- c(ifelse(drawextaxes,3,1),ifelse(drawextaxes,2.3,0.8))
  mtext(c(xlab,ylab),side=c(1,2),line=lab.line,at=c(mean(range(coord[,1])),mean(range(coord[,2]))))
  arrows(coordx1,coordy1,coordx2,coordy2,length=0.06,angle=20,col=col,lwd=lwd)
  if (ident) {
    pos.lab <- function(x1,y1,x2,y2) {
	res <- integer(length(x1))
	for (i in 1:length(x1)) {
	  x2.i <- x2[i]-x1[i]
	  y2.i <- y2[i]-y1[i]
	  res[i] <- if (y2.i<0) {
	    if (x2.i<0) {ifelse(atan(abs(y2.i/x2.i))*180/pi<40,4,3)} else {ifelse(atan(abs(y2.i/x2.i))*180/pi<40,2,3)}
	  } else {
	    if (x2.i<0) {ifelse(atan(abs(y2.i/x2.i))*180/pi<40,4,1)} else {ifelse(atan(abs(y2.i/x2.i))*180/pi<40,2,1)}
	  }
	}
	return(res)
    }
    if (is.null(labels)) {labels <- rownames(coord[as.numeric(pairs)==1,])}
    pos <- pos.lab(coordx1,coordy1,coordx2,coordy2)
    text(coordx1,coordy1,labels,cex=cex,col=col,pos=pos,offset=0.2)
  }
  if (!is.null(main)) {
    main.pos <- match.arg(main.pos)
    xmain <- if (main.pos %in% c("bottomleft","topleft")) {xlim[1]-0.02*diff(xlim)} else {xlim[2]+0.02*diff(xlim)}
    ymain <- if (main.pos %in% c("bottomleft","bottomright")) {ylim[1]} else {ylim[2]}
    adjmain <- if (main.pos %in% c("bottomleft","topleft")) {c(0,NA)} else {c(1,NA)}
    text(xmain,ymain,main,adj=adjmain,cex=main.cex)
  }
  if (legend) {
    if (is.null(legend.lab)) {legend.lab <- "1"}
    if (!is.null(legend.title) && nchar(legend.title)>0) {
	legend(legend.pos,legend.lab,col=legend.col,lty=1,lwd=legend.lwd,bg="white",title=legend.title)
    } else {
	legend(legend.pos,legend.lab,col=legend.col,lty=1,lwd=legend.lwd,bg="white")
    }
  }
  box()
}

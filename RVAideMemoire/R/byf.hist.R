byf.hist <- function(formula,data,sep=FALSE,density=TRUE,xlab=NULL,ylab=NULL) {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call(expand.dots=FALSE)
  m$sep <- m$density <- NULL
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$... <- m$xlab <- m$ylab <- NULL
  mf <- eval(m,parent.frame())
  dname <- c(names(mf)[1],paste(names(mf)[2:ncol(mf)],collapse=":"))
  resp <- mf[,1]
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  if (is.null(xlab)) {xlab <- dname[1]}
  if (is.null(ylab)) {ylab <- "Density"}
  if (sep) {
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))
    par(mfrow=n2mfrow(nlevels(fact)))
    if (density) {
	for (i in 1:nlevels(fact)) {
	  y <- resp[as.numeric(fact)==i]
	  h <- suppressWarnings(hist(y,freq=FALSE,plot=FALSE))
	  plot(0,xlim=range(h$breaks),ylim=c(0,max(h$density)),xlab=xlab,
	    ylab=ylab,main=levels(fact)[i],cex=0)
	  dens <- density(y)
	  col <- col2rgb(palette()[i])
	  col2 <- rgb(col[1,],col[2,],col[3,],alpha=0.4*255,maxColorValue=255)
	  polygon(dens$x,dens$y,col=col2,border=NA)
	  rug(y,col=i)
	}	
    } else {
	for (i in 1:nlevels(fact)) {
	  y <- resp[as.numeric(fact)==i]
	  hist(y,xlab=xlab,main=levels(fact)[i])
	}
    }
  } else {
    if (density) {
	dhist(resp,fac=fact,col=1:nlevels(fact),legend=TRUE,pos.legend="topright",
	  xlab=xlab,ylab=ylab)
    } else {
	nlev <- nlevels(fact)
	couleurs <- 1:nlev
	angles <- seq(20,160,140/(nlev-1))
	dens <- integer(nlev)
	for (i in 1:nlev) {
	  x.i <- resp[as.numeric(fact)==i]
	  dens[i] <- max(hist(x.i,plot=FALSE)$counts)
	}
	x.1 <- resp[as.numeric(fact)==1]
	hist(x.1,xlab=xlab,ylab=ylab,density=10,angle=angles[1],col=couleurs[1],
	  xlim=range(resp),ylim=c(0,max(dens)),main="")
	box()
	for (i in 2:nlev) {
	  x.i <- resp[as.numeric(fact)==i]
	  hist(x.i,xlab="",ylab="",main="",density=10,angle=angles[i],col=couleurs[i],
	    add=TRUE)
	}
	legend("topright",levels(fact),fill=couleurs)
    }
  }
}


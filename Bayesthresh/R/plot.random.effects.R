plot.random.effects <- function(x, by=NULL, interval="confidence", level=0.95, xlab=NULL, ylab=NULL, 
																																xlim=NULL, ylim=NULL, pch=19, 
																																col.points="red", col.seg=gray(0.5), ...)

{
				if(!inherits(x, "random.effects"))
								stop("Use an object of class Bayesthresh")
				if(any(interval=="hpd")){
								if(any(class(x)=="confidence")) stop("HPDinterval is not TRUE in random.effects") 
								if(is.data.frame(x) == TRUE){
												means <- x[,1]
												minimum <- x[,3]
												maximum <- x[,4]
												nomes <- rownames(x)
												dat <- data.frame(means,minimum,maximum,nomes)
												dat <- dat[do.call(order, dat),]
												if(is.null(ylim)) ylim <- c(min(minimum),max(maximum))
												if(is.null(ylab)) ylab = c('Mean')
												plot(dat$means, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, pch=pch, axes=FALSE, 
																	col=col.points, ...)
												axis(1, at=c(1:length(nomes)), labels=dat$nomes, ...)
												dis <- (max(maximum,minimum)-min(minimum,maximum))/10
												axis(2, at=c(round(seq(min(minimum), max(maximum), dis), 2)), ...)
												segments(c(1:length(nomes)),dat$minimum,c(1:length(dat$nomes)),dat$maximum, col=col.seg, ...)
												box()
								}
								else{
												if(is.null(by))
																stop('argument "by" not defined')
												names.x <- names(x)
												pos <- which(names.x == by)
												y <- x[[pos]]
												means <- y[,1]
												minimum <- y[,3]
												maximum <- y[,4]
												nomes <- rownames(y)
												dat <- data.frame(means,minimum,maximum,nomes)
												dat <- dat[do.call(order, dat),]
												if(is.null(ylim)) ylim <- c(min(minimum),max(maximum))
												if(is.null(ylab)) ylab = c('Mean')
												plot(dat$means, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, pch=pch, axes=FALSE, 
																	col=col.points, ...)
												axis(1, at=c(1:length(nomes)), labels=dat$nomes, ...)
												dis <- (max(maximum)-min(minimum))/10
												axis(2, at=c(round(seq(min(minimum), max(maximum), dis), 2)), ...)
												segments(c(1:length(nomes)),dat$minimum,c(1:length(dat$nomes)),dat$maximum, col=col.seg, ...)
												box()

								}
				}

				else {
								if(is.data.frame(x) == TRUE)
								{
												names.x <- names(x)
												means <- x[,1]
												quant <-  abs(qnorm((1-level)/2,mean=x[,1], sd=x[,2]))
												minimum <- means-(x[,2]*quant)
												maximum <- means+(x[,2]*quant)
												nomes <- rownames(x)
												dat <- data.frame(means,minimum,maximum,nomes)
												dat <- dat[do.call(order, dat),]
												if(is.null(ylim)) ylim <- c(min(minimum),max(maximum))
												if(is.null(ylab)) ylab = c('Mean')
												plot(dat$means, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, pch=pch, axes=FALSE, 
																	col=col.points, ...)
												axis(1, at=c(1:length(nomes)), labels=dat$nomes, ...)
												dis <- (max(maximum,minimum)-min(minimum,maximum))/10
												axis(2, at=c(round(seq(min(minimum), max(maximum), dis), 2)), ...)
												segments(c(1:length(nomes)),dat$minimum,c(1:length(dat$nomes)),dat$maximum, col=col.seg, ...)
												box()
								}
								else {
												if(is.null(by))
																stop('argument "by" not defined')
												names.x <- names(x)
												pos <- which(names.x == by)
												y <- x[[pos]]
												means <- y[,1]
												quant <-  abs(qnorm((1-level)/2,mean=y[,1], sd=y[,2]))
												minimum <- means-(y[,2]*quant)
												maximum <- means+(y[,2]*quant)
												nomes <- rownames(y)
												dat <- data.frame(means,minimum,maximum,nomes)
												dat <- dat[do.call(order, dat),]
												if(is.null(ylim)) ylim <- c(min(minimum),max(maximum))
												if(is.null(ylab)) ylab = c('Mean')
												plot(dat$means, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, pch=pch, axes=FALSE, 
																	col=col.points, ...)
												axis(1, at=c(1:length(nomes)), labels=dat$nomes, ...)
												dis <- (max(maximum)-min(minimum))/10
												axis(2, at=c(round(seq(min(minimum), max(maximum), dis), 2)), ...)
												segments(c(1:length(nomes)),dat$minimum,c(1:length(dat$nomes)),dat$maximum, col=col.seg, ...)
												box()
								}
				}
}


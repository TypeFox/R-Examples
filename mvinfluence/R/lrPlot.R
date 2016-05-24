## LR influence plot: plot of log studentized residual^2 vs log(h/p*(1-h))
## such that contours of equal Cook's distance fall on diagonal lines
## with slope = -1
# Following McCulloch & Meeter, Technometrics, 1983, 25, 152-155
# modified from car::influencePlot

lrPlot <- function(model, ...){
    UseMethod("lrPlot")
    }

lrPlot.lm <- function(model, scale=12,  
		xlab="log Leverage factor [log H/p*(1-H)]", ylab="log (Studentized Residual^2)",
		xlim=NULL, ylim, 
    labels, 
    id.method = "noteworthy",
    id.n = if(id.method[1]=="identify") Inf else 0, 
    id.cex=1, id.col=palette()[1],
    ref=c("h", "v", "d", "c"),                     # reference lines
    ref.col="gray", ref.lty=2, ref.lab=TRUE,       # args for reference lines
     ...){ 

	hatval <- hatvalues(model)
	rstud <- rstudent(model)
	p <- length(coef(model))

	Hfun <- function(hat) log(hat/(p*(1-hat)))
	
	L <- Hfun(hatval)
	R <- log(rstud^2)
	if (missing(labels)) labels <- names(rstud)
	cook <- sqrt(CookD<-cooks.distance(model))
	scale <- scale/max(cook, na.rm=TRUE)
	n <- sum(!is.na(rstud))
#	cutoff <- sqrt(4/(n - p))

	# allow 6 log steps by default to avoid small values swamping the plot
	if (missing(ylim)) ylim <- rev(c(ymax<-ceiling(max(R)), ymax-6))

	plot(L, R, xlab=xlab, ylab=ylab, type='n', xlim=xlim, ylim=ylim,  ...)
	points(L, R, cex=pmax(.05, scale*cook), ...)

	nomit <- sum(R < ylim[1])
	if (nomit>0) cat("Note:",nomit, "points less than R=", ylim[1], "have been clipped to preserve resolution\n")

	# diagonal lines of slope -1, indicating constant leverage
	xmin <- floor(par("usr")[1])
	xmax <- floor(par("usr")[2])
	# FIXME: bit of a kludge in calculating intercepts of diagonal lines, but does no harm
	if ("d" %in% ref) for(a in (2*xmin):xmax) abline(a=a, b=-1, col=ref.col, lty=ref.lty)
	if ("c" %in% ref) {
		cref <- c(.95, .99)
		for(i in seq_along(cref)) {
			a <- -log(qf(cref[i],p,n-p))
			abline(a=a, b=-1, lwd=2, col="red")
			if(ref.lab) text(xmax, a-xmax, cref[i], adj=c(1,1), cex=0.9, col="red")
			}
		}

	# high leverage points, on this scale
	if ("v" %in% ref) {
		hL <- Hfun(c(2, 3)*p/n)
		ymin <- par("usr")[3]
		abline(v=hL, lty=ref.lty, col=ref.col)
		if(ref.lab) text(hL, ymin, c("2p/n", "3p/n"), pos=3, cex=0.9)
		}
	# high rstud^2
	if ("h" %in% ref) {
		abline(h=log(qf(.95, 1, n-p)), lty=ref.lty, col=ref.col)
		if(ref.lab) text(par("usr")[1], log(qf(.95, 1, n-p)), "0.95", adj=c(0,0), cex=0.9)
		}
	
	if(id.method != "identify"){
	   which.rstud <- order(abs(rstud), decreasing=TRUE)[1:id.n]
	   which.cook <- order(cook, decreasing=TRUE)[1:id.n]
	   which.hatval <- order(hatval, decreasing=TRUE)[1:id.n]
	   which.all <- union(which.rstud, union(which.cook, which.hatval))
	   id.method <- which.all
	   }

	noteworthy <- car:::showLabels(L, R, labels=labels, id.method=id.method, 
    id.n=id.n, id.cex=id.cex, id.col=id.col)
  	if (length(noteworthy > 0))
	res <- data.frame(Rstudent=rstud, Hat=hatval, CookD=CookD, L, R)[noteworthy,]
	return(res)
  }


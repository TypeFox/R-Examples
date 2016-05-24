influencePlot.mlm <-
function(model, scale=12, type=c("cookd", "stres", "LR"),
		infl=mlm.influence(model, do.coef = FALSE), FUN=det,
    fill=TRUE, fill.col="blue",
    fill.alpha.max=0.5,
    labels, 
    id.method = "noteworthy",
    id.n = if(id.method[1]=="identify") Inf else 0, 
    id.cex=1, id.col=palette()[1],
    ref.col="gray", ref.lty=2, ref.lab=TRUE,       # args for reference lines
    ...)
{
	
	m <- infl$m
	df <- as.data.frame(infl, FUN=FUN, funnames=FALSE)
	CookD <- df$CookD
	H <- df$H
	Q <- df$Q
	L <- df$L
	R <- df$R
	if(missing(labels)) labels <- rownames(df)
	
	p <- nrow(coef(model))
	n <- length(infl$H)
	scale=scale/max(sqrt(CookD))
	if (fill) cols <- heplots:::trans.colors(fill.col, alpha=fill.alpha.max * CookD/(max(CookD)))

	if(id.method != "identify"){
	  which.rstud <- order(R, decreasing=TRUE)[1:id.n]
		which.cook <- order(CookD, decreasing=TRUE)[1:id.n]
		which.hat <- order(H, decreasing=TRUE)[1:id.n]
		id.method <- inf <- Reduce(union, list(which.rstud,which.cook, which.hat))
	}
	
	type <- match.arg(type)
	if (type=='cookd') {
		plot(H, CookD, xlab="Hat value", ylab="Cook's D", cex=scale*CookD, ...)
		if (fill) points(H, CookD, cex=scale*CookD, pch=16, col=cols)
#		abline(h=qf(.95, 1, n-p), lty=ref.lty, col=ref.col)
		abline(h=4/(n-p), lty=ref.lty, col=ref.col)
		abline(v=c(2, 3)*p/n, lty=ref.lty, col=ref.col)
		noteworthy <- car:::showLabels(H, CookD, labels=labels, id.method=id.method, 
    	id.n=id.n, id.cex=id.cex, id.col=id.col)
	}
	else if (type=='stres') {
		plot(H, Q, xlab="Hat value", ylab="Squared Residual", cex=scale*CookD, ...)
		if (fill) points(H, Q, cex=scale*CookD, pch=16, col=cols)
		abline(v=c(2, 3)*p/n, lty=ref.lty, col=ref.col)
		noteworthy <- car:::showLabels(H, Q, labels=labels, id.method=id.method, 
    	id.n=id.n, id.cex=id.cex, id.col=id.col)
	}
	else if (type=='LR') {
		logL <- log(L)
		logR <- log(R)
		plot(logL, logR, xlab="log Leverage", ylab="log Residual", cex=scale*CookD, ...)
		if (fill) points(logL, logR, cex=scale*CookD, pch=16, col=cols)
		xmin <- floor(par("usr")[1])
		xmax <- ceiling(par("usr")[2])
		# FIXME: bit of a kludge in calculating intercepts of diagonal lines
		for(a in (2*xmin):xmax) abline(a=a, b=-1, col=ref.col, lty=ref.lty)	
		noteworthy <- car:::showLabels(logL, logR, labels=labels, id.method=id.method, 
    	id.n=id.n, id.cex=id.cex, id.col=id.col)
		}
#browser()

	  if (length(noteworthy > 0)) {
	  	## FIXME:  m>1 car::showLabels returns the labels, not the indices
	  	if (m>1) noteworthy <- inf
	  	res <- data.frame(H, Q, CookD, L, R)[noteworthy,]
			return(res)
		}

}

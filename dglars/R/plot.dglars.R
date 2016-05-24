plot.dglars <- function(x,k=c("BIC","AIC"),complexity=c("df","gdf"),g.gof=NULL,...){
	if(is.numeric(k)){
		if(k<=0) stop("k must be greater than zero")
		knm <- "GoF"
	}
	else{
		knm <- match.arg(k)
		k <- ifelse(knm == "BIC",log(dim(x$X)[1]),2)
	}
	complexity <- match.arg(complexity)
	if(x$family=="poisson" & complexity == "gdf"){
		complexity <- "df"
		warning("'complexity' was set equal to 'df' because for Poisson regression the\n model complexity can be approximated by the number of nonzero coefficients")
	}
	n <- dim(x$X)[1]
	beta <- x$beta
	dev <- x$dev
	g <- x$g
	g.action <- g[x$action!=""]
	if(is.null(g.gof)){
		if(complexity == "df")	df <- x$df
		else df <- gdf(x)
		gof <- dev + k * df
		g.gof <- g[which.min(gof)]
		plot(g,gof,xlab=expression(gamma),ylab=knm,type="n",main="Model Selection Criterion")
		abline(v=g.action,lty=2,col=8)
		abline(v=g.gof,lty=2,col=2)
		points(g,gof,xlab=expression(gamma),ylab=knm,type="o",pch=20,lty=2,...)
		op <- par(ask=dev.interactive())
		
	} else knm <- "g.gof"
	matplot(g,t(beta[-1,]),col=1,type="n",xlab=expression(gamma),ylab="Regression Coefficients",main="Coefficients Path")
	abline(v=g.action,lty=2,col=8)
	abline(v=g.gof,lty=2,col=2)	
	matpoints(g,t(beta[-1,]),col=1,type="l",lty=1,...)
	axis(3,g.gof,knm,padj=1)
	if(!is.null(g.gof)) op <- par(ask=dev.interactive())
	if(x$control$algorithm=="pc"){
		ru <- x$ru
		matplot(g,t(ru),col=1,type="n",xlab=expression(gamma),ylab="Rao Score Statistics",main="Rao Score Path")
		abline(v=g.action,lty=2,col=8)
		abline(v=g.gof,lty=2,col=2)
		matpoints(g,t(ru),col=1,type="l",lty=1,...)
		if(length(dev)!=1) axis(3,g.gof,knm,padj=1)
	}
	par(op)
}
"nlsConfRegions"<-function(nls, length=1000, exp=1.5){
	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	data	<- eval(nls$call$data, sys.frame(0))
	np		<- length(coef(nls))
	nl		<- nrow(data)
	vardep <- all.vars(formula(nls)[[2]])
	varindep <- intersect(all.vars(formula(nls)[[3]]), colnames(data))
	"formula2function"<-function(formu){
		arg1		<- all.vars(formu)
		arg2		<- vector("list",length(arg1))
		names(arg2)	<- arg1
		Args		<- do.call("alist",arg2)
		fmodele		<- as.function(c(Args,formu))
		return(fmodele)
	}

	fmodele		<- formula2function(formula(nls)[[3]])
	scer <- sum(residuals(nls)^2)
	scer95 <- scer * (1 + (np/(nl-np)) * qf(p=0.95, df1=np, df2=(nl-np), lower.tail=TRUE))
	student95	<- qt(0.975, df=nl-np)
	bornes <- cbind(coef(nls)-student95*exp*(summary(nls)$parameters)[,2], coef(nls)+student95*exp*(summary(nls)$parameters)[,2])
	colnames(bornes) <- c("Lower", "Upper")
	tirage	<- vector(length=np)
	names(tirage)	<- row.names(summary(nls)$parameters)
	tab	<- matrix(ncol=np, nrow=0)
	rss	<- vector(length=0)

	cat("  ")
	while(nrow(tab)<length){
		tirage <- apply(bornes, 1, function(z) runif(n=1, min=z[1], max=z[2]))
		listparavar	<- c(tirage, data[varindep])
		predict	<- do.call("fmodele", listparavar)
		rss1	<- sum((predict-data[,vardep])^2)
		if(rss1 < scer95){
			tenth	<- floor(100*nrow(tab)/length)
			tab	<- rbind(tab,tirage)
			rss	<- c(rss,rss1)
			if(tenth!=floor(100*nrow(tab)/length)){
				if(tenth<11){cat("\b\b",tenth,"%",sep="")}
				else{cat("\b\b\b",tenth,"%",sep="")}
			}
		}
    }
    rownames(tab) <- 1:nrow(tab)
	cat("\b\b\b100%\a")
	cat("\n Confidence regions array returned \n")
	listcr <- list(cr=tab, rss=rss, rss95=scer95, bounds=bornes)
	class(listcr) <- "nlsConfRegions"
	return(listcr)
}


"plot.nlsConfRegions"<-function(x, bounds=FALSE, ask=FALSE, ...){
	if (!inherits(x, "nlsConfRegions"))
		stop("Use only with 'nlsConfRegions' objects")
	np <- ncol(x$cr)
	def.par <- par(no.readonly = TRUE)
	if(ask) par(ask=TRUE,mar=c(4,4,3,1))
	if(!ask){
		lay <- lower.tri(matrix(0,(np-1),(np-1)), TRUE)
		lay[which(lay, TRUE)] <- 1:choose(np,2)
		layout(lay)
		par(mar=c(5,4,0.2,0.2))
	}
	for(i in 1:(np-1))
		for(j in (i+1):np){
			if(!bounds) plot(x$cr[,i], x$cr[,j], pch="+", xlab=colnames(x$cr)[i], ylab=colnames(x$cr)[j])
			else{
				xrange <- range(c(x$cr[,i], x$bounds[i,]))
				yrange <- range(c(x$cr[,j], x$bounds[j,]))
				plot(x$cr[,i], x$cr[,j], pch="+", xlab=colnames(x$cr)[i], ylab=colnames(x$cr)[j], xlim=xrange, ylim=yrange)
				abline(h=x$bounds[j,], v=x$bounds[i,], col="red", lty=2)
			}
		}

	par(def.par)	
}

"print.nlsConfRegions" <- function (x, ...) {
	if (!inherits(x, "nlsConfRegions"))
		stop("Use only with 'nlsConfRegions' objects")
	cat("Beale's 95 percent confidence regions\n")
	cat("\n")
	sumry <- array("", c(2, 4), list(1:2, c("vector", "length", "mode", "content")))
	sumry[1, ] <- c("$rss", length(x$rss), mode(x$rss), "Residual sums of squares")
	sumry[2, ] <- c("$rss95", length(x$rss95), mode(x$rss95), "95 percent RSS threshold")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(2, 4), list(1:2, c("data.frame", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$cr", nrow(x$cr), ncol(x$cr), "Sample drawn in the confidence region")
	sumry[2, ] <- c("$bounds", nrow(x$bounds), ncol(x$bounds), "Bounds of the drawing hypercube")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}

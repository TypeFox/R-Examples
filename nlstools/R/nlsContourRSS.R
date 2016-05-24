"nlsContourRSS"<-function(nls, lseq=100, exp=2){
	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	"formula2function"<-function(formu){
		arg1		<- all.vars(formu)
		arg2		<- vector("list",length(arg1))
		names(arg2)	<- arg1
		Args		<- do.call("alist",arg2)
		fmodele		<- as.function(c(Args,formu))
		return(fmodele)
	}

	"sce" <- function(para1, para2, i, j, para=lestimates, vari=data[varindep], resp=data[,vardep]){
		para[[i]] <- para1; para[[j]] <- para2
		lvari <- as.list(vari)
		paraVar <- c(para, lvari)
		pred <- do.call("fmodele", paraVar)
		sum((resp-pred)^2)
	}

	vSCE <- Vectorize(sce, c("para1", "para2"))

	data	<- eval(nls$call$data, sys.frame(0))
	np		<- length(coef(nls))
	ncomb	<- choose(np, 2)
	vardep <- all.vars(formula(nls)[[2]])	
	varindep <- intersect(all.vars(formula(nls)[[3]]), colnames(data))
	nl		<- nrow(data)
	fmodele		<- formula2function(formula(nls)[[3]])
	estimates <- coef(nls)
	lestimates <- as.list(estimates)
	scer <- sum(residuals(nls)^2)
	scer95		<- scer * (1 + (np/(nl-np)) * qf(p=0.95, df1=np, df2=(nl-np), lower.tail=TRUE))
	student95	<- qt(0.975, df=nl-np)
	bornes <- cbind(coef(nls)-student95*exp*(summary(nls)$parameters)[,2], coef(nls)+student95*exp*(summary(nls)$parameters)[,2])
	rownames(bornes) <- names(estimates)
	lsce <- list()

	seqPara <- apply(bornes, 1, function(z) seq(from=z[1], to=z[2], length=lseq))
	esti2 <- as.list(estimates)

	for(i in 1:(np-1)){
		for(j in (i+1):np){
			tenth <- floor(100*length(lsce)/choose(np,2))
			m12 <- outer(seqPara[,i], seqPara[,j], vSCE, i, j)
			lsce <- c(lsce, list(m12))
			if(tenth!=floor(100 * length(lsce) / choose(np,2))){
				if(tenth<11){cat("\b\b",tenth,"%",sep="")}
				else{cat("\b\b\b",tenth,"%",sep="")}
			}
		}
	}
	cat("\b\b\b100%\a")
	cat("\n RSS contour surface array returned \n")
	llsce <- lapply(lsce, log)
	listsce <- list(seqPara=seqPara, lrss=llsce, lrss95=log(scer95))
	class(listsce) <- "nlsContourRSS"
	return(listsce)
}


"plot.nlsContourRSS"<-function(x, nlev=0, col=TRUE, col.pal=terrain.colors(100), ask=FALSE, useRaster=TRUE, ...){
	if (!inherits(x, "nlsContourRSS"))
		stop("Use only with 'nlsContourRSS' objects")

	np <- ncol(x$seqPara)
	paranames <- colnames(x$seqPara)
	count <- 0
 	def.par <- par(no.readonly = TRUE)
 	
 	if(ask) par(ask=TRUE)
	if(!ask){
		lay <- lower.tri(matrix(0,(np-1),(np-1)), TRUE)
		lay[which(lay, TRUE)] <- 1:choose(np,2)
		layout(lay)
		par(mar=c(5,4,0.2,0.2))
	}
	for(i in 1:(np-1))
		for(j in (i+1):np){
			count <- count + 1
			if(col){
				image(x$seqPara[,i], x$seqPara[,j], x$lrss[[count]], xlab=paranames[i], ylab=paranames[j], col=col.pal, useRaster=useRaster)
				if(nlev>0) contour(x$seqPara[,i], x$seqPara[,j], x$lrss[[count]], labels="", nlevels=nlev, add=TRUE)
				contour(x$seqPara[,i], x$seqPara[,j], x$lrss[[count]], labels="", levels=x$lrss95, lty=3, col="red", lwd=2, add=TRUE)
			}
			else{
				if(nlev>0){
					contour(x$seqPara[,i], x$seqPara[,j], x$lrss[[count]], labels="", nlevels=nlev, xlab=paranames[i], ylab=paranames[j])
					contour(x$seqPara[,i], x$seqPara[,j], x$lrss[[count]], labels="", levels=x$lrss95, lty=3, col="red", lwd=2, add=TRUE)
				}
				else contour(x$seqPara[,i], x$seqPara[,j], x$lrss[[count]], labels="", levels=x$lrss95, lty=3, col="red", lwd=2)
			}
		}
	par(def.par)	
}

"print.nlsContourRSS" <- function (x, ...) {
	if (!inherits(x, "nlsContourRSS"))
		stop("Use only with 'nlsContourRSS' objects")
	cat("RSS surface contour\n")
	cat("\n")
	cat("$lrss95 (95 percent RSS threshold): ", signif(x$lrss95,3),"\n")
	cat("\n")
	sumry <- array("", c(1, 4), list(1:1, c("matrix", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$seqPara", nrow(x$seqPara), ncol(x$seqPara), "Sequence of parameters")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(1, 5), list(1:1, c("list", "length", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$lrss", length(x$lrss), nrow(x$lrss[[1]]), ncol(x$lrss[[1]]), "RSS")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}

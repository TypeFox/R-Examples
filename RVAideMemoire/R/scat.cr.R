scat.cr <- function(dudi.obj,axis=1) {
  if (!"acm"%in%class(dudi.obj) & !"mix"%in%class(dudi.obj)) {stop("unknown analysis")}
  cr <- dudi.obj$cr
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  par(mar=c(7,4,4,1))
  if ("acm" %in% class(dudi.obj)) {
    barplot(cr[,axis],names.arg=rownames(cr),las=2,ylab=paste("Correlation ratio to axis",axis))
  } else if ("mix" %in% class(dudi.obj)) {
    index <- dudi.obj$index
    col <- rep("",length(index))
    col[which(index=="q")] <- "salmon"
    col[which(index=="f")] <- "lightblue"
    if ("o" %in% levels(index)) {
	col[which(index=="o")] <- "mediumpurple"
    }
    barplot(cr[,axis],names.arg=rownames(cr),las=2,col=col,ylab=paste("Value relatively to axis",axis))
    if (nlevels(index)==3) {
	legend("topleft",c("Squared corr coef","Squared multiple corr coef","Correlation ratio"),
	  fill=c("salmon","mediumpurple","lightblue"),cex=0.8)
    } else if (nlevels(index)==2) {
	if("q" %in% levels(index) & "o" %in% levels(index)) {
	  legend("topleft",c("Squared corr coef","Squared multiple corr coef"),
	    fill=c("salmon","mediumpurple"),cex=0.8)
	} else if ("q" %in% levels(index) & "f" %in% levels(index)) {
	  legend("topleft",c("Squared corr coef","Correlation ratio"),
	    fill=c("salmon","lightblue"),cex=0.8)
	} else if ("o" %in% levels(index) & "f" %in% levels(index)) {
	  legend("topleft",c("Squared multiple corr coef","Correlation ratio"),
	    fill=c("mediumpurple","lightblue"),cex=0.8)
	}
    } else {
	if(levels(index)=="q") {
	  legend("topleft","Squared corr coef",fill="salmon",cex=0.8)
	} else if (levels(index)=="o") {
	  legend("topleft","Squared multiple corr coef",fill="mediumpurple",cex=0.8)
	} else if (levels(index)=="f") {
	  legend("topleft","Correlation ratio",fill="lightblue",cex=0.8)
	}
    }
  }
}

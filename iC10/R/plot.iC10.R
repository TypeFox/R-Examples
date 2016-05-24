plot.iC10 <-
function(x, sample.name=1, newdata=NULL,...) {
    par(no.readonly=TRUE)
    oldpar <- par(no.readonly=TRUE)
    if (!is.numeric(sample.name)) sample.name <- which(names(x$class) == sample.name)
    coliCluster <- c('#FF5500', '#00EE76', '#CD3278', '#00C5CD', '#8B0000',
                     '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
    cn.features <- NULL
    exp.features <- NULL
    if (attr(x, "classifier.type") != "Exp") {
          cn.features <- x$map.cn
          cn.features$Type <- "CN"
      }
      if (attr(x, "classifier.type") != "CN") {
          exp.features <- x$map.exp[,-c(9,10)]
          exp.features$Type <- "Exp"
   }
    max.grey.square <- nrow(cn.features)
    features <- rbind(cn.features, exp.features)
    if (attr(x, "ref") == "hg18") {
        features$CHROM <- features$chromosome_name_hg18
    }
    if (attr(x, "ref") == "hg19") {
        features$CHROM <- features$chromosome_name_hg19
    }
    Pos <- c(1, which(diff(features$CHROM)!=0), nrow(features))
    text.pos <- Pos + c(diff(Pos)/2, 0)
    text.pos <- text.pos[-length(text.pos)]
    if (!is.null(newdata)) {
        tmp <- rbind(newdata$CN, newdata$Exp)
        for (i in 1:10) {
            if (length(which(x$class==i))>0) {
                x$centroids[,i] <- apply(tmp[,which(x$class==i),drop=F], 1, mean)
            }
        }
    }
    par(mfrow=c(6, 2), mar=c(2.5, 4, 2, 2), oma=c(2, 2, 2, 2), no.readonly=TRUE)
    for (i in 1:10) {
        plot(x$centroids[,i], type="h", col=coliCluster[i], xlab="", ylab=paste("iC", i), axes=F, 
	ylim=range(x$centroids))
        axis(2)
        bottom.axis.pos <- seq(from=1, by=2, length=length(text.pos))
        top.axis.pos <- seq(from=2, by=2, length=length(text.pos))
        axis(1, at=text.pos[bottom.axis.pos], features$CHROM[Pos[-1]][bottom.axis.pos])
        axis(3, at=text.pos[top.axis.pos], features$CHROM[Pos[-1]][top.axis.pos])
        box()
        if (attr(x, "classifier.type") != "Exp") {
            polygon(x=rep(c(0, max.grey.square), c(2, 2)),
                    y=c(range(x$centroids[,i]),rev(range(x$centroids[,i]))),
                    density=NA, border="lightgrey", col="lightgrey", ylab=paste("iC", i))
        }
        points(x$centroids[,i], type="h", col=coliCluster[i], xlab="")
        abline(v=Pos, lty=2)
    }
    if (is.null(newdata)) {
        plot(x$fitted[,sample.name], type="h", col="black", xlab="", ylab="Prediction", axes=F)
        axis(2)
        axis(1, at=text.pos[bottom.axis.pos], features$CHROM[Pos[-1]][bottom.axis.pos])
        axis(3, at=text.pos[top.axis.pos], features$CHROM[Pos[-1]][top.axis.pos])
        if (attr(x, "classifier.type") != "Exp") {
            polygon(x=rep(c(0, max.grey.square), c(2, 2)),
                    y=c(range(x$fitted[,sample.name]),rev(range(x$fitted[,sample.name]))),
                    density=NA, border="lightgrey", col="lightgrey")
        }
        points(x$fitted[,sample.name], type="h", col="black")
        abline(v=Pos, lty=2)
        barplot(x$posterior[sample.name,], col=coliCluster, ylim=c(0,1), ylab="Prob.classif.")
    } else {
        plot(0,0, type="n", axes=F, xlab="", ylab="")
        barplot(table(x$class), col=coliCluster, ylab="Frequency")
        }
    par(oldpar, no.readonly=TRUE)
}

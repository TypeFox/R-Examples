compare <- function(obj, iC10=1:10, newdata, name.test="Test",...) {
	UseMethod("compare")
}
compare.iC10 <-
function(obj, iC10=1:10, newdata, name.test="Test",...) {
    par(no.readonly=TRUE)
    oldpar <- par(no.readonly=TRUE)
    coliCluster <- c('#FF5500', '#00EE76', '#CD3278', '#00C5CD', '#8B0000',
                     '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
    cn.features <- NULL
    exp.features <- NULL
    if (attr(obj, "classifier.type") != "Exp") {
        cn.features <- newdata$map.cn
	if ("Synonyms_0" %in% colnames(cn.features)) {
	   cn.features <- cn.features[,-which(colnames(cn.features) %in% c("Synonyms_0", "Gene.Chosen"))]
         }
        cn.features$Type <- "CN"
	features <- cn.features
    }
    if (attr(obj, "classifier.type") != "CN") {
        exp.features <- newdata$map.exp
        exp.features$Type <- "Exp"
	features <- exp.features
    }
    max.grey.square <- nrow(cn.features)
    if (!is.null(cn.features) & !is.null(exp.features)) {
        common.cols <- intersect(colnames(cn.features), colnames(exp.features))
       features <- rbind(cn.features[,common.cols], exp.features[,common.cols])
    } 

    if (attr(obj, "ref") == "hg18") {
        features$CHROM <- features$chromosome_name_hg18
    }
    if (attr(obj, "ref") == "hg19") {
        features$CHROM <- features$chromosome_name_hg19
    }
    Pos <- c(1, which(diff(features$CHROM)!=0), nrow(features))
    text.pos <- Pos + c(diff(Pos)/2, 0)
    text.pos <- text.pos[-length(text.pos)]
    all.data <- rbind(newdata$CN, newdata$Exp)
    par(mfrow=c(length(iC10), 2), mar=c(2.5, 4, 5, 2), oma=c(2, 2, 2, 2), no.readonly=TRUE)
    for (i in iC10) {
        plot(obj$centroids[,i], type="h", col=coliCluster[i], xlab="", ylab=paste("iC", i), axes=F,
             main=paste("Training data. iC", i, sep=""),...)
        axis(2)
        bottom.axis.pos <- seq(from=1, by=2, length=length(text.pos))
        top.axis.pos <- seq(from=2, by=2, length=length(text.pos))
        axis(1, at=text.pos[bottom.axis.pos], features$CHROM[Pos[-1]][bottom.axis.pos])
        axis(3, at=text.pos[top.axis.pos], features$CHROM[Pos[-1]][top.axis.pos])
        box()
        if (attr(obj, "classifier.type") != "Exp") {
            polygon(x=rep(c(0, max.grey.square), c(2, 2)),
                    y=c(range(obj$centroids[,i]),rev(range(obj$centroids[,i]))),
                    density=NA, border="lightgrey", col="lightgrey", ylab=paste("iC", i))
        }
        points(obj$centroids[,i], type="h", col=coliCluster[i], xlab="")
        abline(v=Pos, lty=2)
        if (length(which(obj$class==i))>0) {
            test <- unclass(apply(all.data[,which(obj$class==i),drop=F], 1, mean))
            plot(test, type="h", col=coliCluster[i], xlab="", ylab=paste("iC", i), axes=F,
                 main=paste(name.test, " data. iC", i, sep=""),...)
            axis(2)
            bottom.axis.pos <- seq(from=1, by=2, length=length(text.pos))
            top.axis.pos <- seq(from=2, by=2, length=length(text.pos))
            axis(1, at=text.pos[bottom.axis.pos], features$CHROM[Pos[-1]][bottom.axis.pos])
            axis(3, at=text.pos[top.axis.pos], features$CHROM[Pos[-1]][top.axis.pos])
            box()
            if (attr(obj, "classifier.type") != "Exp") {
                polygon(x=rep(c(0, max.grey.square), c(2, 2)),
                        y=c(range(obj$centroids[,i]),rev(range(obj$centroids[,i]))),
                        density=NA, border="lightgrey", col="lightgrey", ylab=paste("iC", i))
            }
            points(test, type="h", col=coliCluster[i], xlab="")
            abline(v=Pos, lty=2)
        } else {
	plot(0,0, type="n", xlab="", ylab="", axes=F)
}
    }
    par(oldpar, no.readonly=TRUE)
}

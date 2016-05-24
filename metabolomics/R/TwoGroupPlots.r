TwoGroupPlots <- function(datamat, tstats, foldchanges, pvalues, 
    padjmethod="BH", fcutoff=log(2), pcutoff=0.05, cexval=0.7)
{
    # Sanity checks
    datamat<-as.matrix(datamat)
    if (mode(datamat)!="numeric" ) {
        stop("datamat must be a numerical data matrix")
    }
    if (class(tstats)!="numeric") {
        stop("tstats must be a vector")
    }
    if (class(foldchanges)!="numeric") {
        stop("foldchanges must be a vector")
    }
    if (class(pvalues)!="numeric") {
        stop("pvalues must be a vector")
    }
    
    avg <- rowMeans(t(datamat), na.rm=TRUE)
    names(tstats)<-names(foldchanges)<-names(pvalues)<-names(avg)
    
    padj <- p.adjust(pvalues, padjmethod)
    cutoff <- which(padj < pcutoff & abs(foldchanges) > fcutoff)
    mark_t <- tstats[cutoff]
    mark_coeff <- foldchanges[cutoff]
    mark_avg <- avg[cutoff]
    signames_up <- names(avg)[which(padj < pcutoff & foldchanges > fcutoff)]
    signames_down <- names(avg)[which(padj < pcutoff & foldchanges < -fcutoff)]
    signames_all <- names(avg)[cutoff]
    
    op <- par(mfrow=c(1, 2))
    #t-statistic vs average plot
    plot(avg, tstats, xlab="Average log abundance", ylab="t statistic",
        col="#3366aa", lwd=2
    )
    points(mark_avg, mark_t, col="#ee3333", pch=16)
    abline(h=0, lty=2, col="#3366aa")
    
    # Labels for average log abundance plot
    for (met in signames_all) {
        xval <- avg[[met]]
        yval <- tstats[[met]]
        xrange <- range(as.numeric(avg), na.rm=TRUE)
        yrange <- range(as.numeric(tstats), na.rm=TRUE)
        # Left plot
        if (yval < (0.95 * max(yrange))) {
            if (xval < (0.75 * max(xrange))) {
                posval <- 4          # right if it's neither
            } else {
                posval <- 2          # left if > 0.75 max(x values)
            }
        } else {
            posval <- 1              # bottom if > 0.95 max(y values)
        }
        text(xval, yval, met, #pos=4,
             pos=posval,
             cex=cexval
        )
    }
    
    #se versus fold change
    plot((foldchanges/tstats), foldchanges, ylab="Fold changes", 
        xlab="Standard errors", col="#3366aa", lwd=2
    )
    points(mark_coeff/mark_t, mark_coeff, col="#ee3333", pch=16)
    abline(h=0, lty=2, col="#3366aa")
    
    for (met in signames_all) {
        xval <- foldchanges[[met]] / tstats[[met]]
        yval <- foldchanges[[met]]
        xrange <- range(as.numeric(foldchanges/tstats), na.rm=TRUE)
        yrange <- range(as.numeric(foldchanges), na.rm=TRUE)
        # Right plot
        if (yval < (0.95 * max(yrange))) {
            if (xval < (0.75 * max(xrange))) {
                posval <- 4          # right if it's neither
            } else {
                posval <- 2          # left if > 0.75 max(x values)
            }
        } else {
            posval <- 1              # bottom if > 0.95 max(y values)
        }
        text(xval, yval, met, #pos=4, 
            pos=posval,
            cex=cexval
        )
    }
    # Reset graphic parameters
    par(op)
    return(
        list(IncreasedMets=signames_up, 
            DecreasedMets=signames_down, 
            DifferentialMets=signames_all
        )
    )
}

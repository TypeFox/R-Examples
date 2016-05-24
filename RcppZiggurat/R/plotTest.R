
plotAll <- function(data) {
    op <- par(mfrow=c(2,3), mar=c(3,3,3,1), oma=c(1,0,2,0))
    for (i in seq_len(ncol(data))) {
        plotTest(data[,i], colnames(data)[i])
    }

    title(paste(attr(data, "testtype"), "test results"), line=0, outer=TRUE, cex.main=2)
    txt <- paste0("Draws:",       attr(data, "draws"),
                  " Repeats: ",   attr(data, "repeats"),
                  " Seed: ",      attr(data, "seed"),
                  " Created at: ",attr(data, "created"),
                  " Version: ",   attr(data, "version"))
    mtext(text=txt, side=1, outer=TRUE, line=-0.5,
          adj=0.02, las=1, cex=0.66)

    par(op)
    invisible(NULL)
}

plotTest <- function(v, g) {
    pks <- ks.test(v, "punif", 0, 1, exact=TRUE)$p.value
    pw <- wilcox.test(v, mu=0.5)$p.value

    plot(ecdf(v), verticals=TRUE, do.p=FALSE,
         xlim=c(0,1), ylim=c(0,1), main=g)
    segments(0,0,1,1, col='darkgray', lty="dotted")
    legvec <- c(paste("pKS:", round(pks, digits=4)),
                paste("pWil.:", round(pw, digits=4)))
    legend(x=-0.1, y=1.0, lty=NULL, bty="n", legend=legvec)
    invisible(NULL)
}

plotChiSq <- function(res, verbose=FALSE) {
    bins <- attr(res, "bins")
    cval <- qchisq(0.95, bins-1)
    if (verbose) {
        cat(sprintf("Critical one-sided 95%% value is %f\n", cval))
        cat(sprintf("Actual chisq(%d) values\n", bins))
        print(tail(res,1))
    }

    op <- par(mfrow=c(2,3), mar=c(3,3,3,1), oma=c(1,0,2,0))
    k <- ncol(res)
    yrange <- c(0,max(max(res[,-1]), cval))
    for (i in 2:k) {
        plot(res[,1], res[,i], type='l', ylim=yrange,
             main=colnames(res)[i])
        abline(h=cval, col="darkgrey", lty="dotted")
    }
    title("Chi-square test results", line=0, outer=TRUE, cex.main=2)
    txt <- paste0("Total draws:", attr(res, "draws"),
                  " Bins: ",      attr(res, "bins"),
                  " Seed: ",      attr(res, "seed"),
                  " Steps: ",     attr(res, "steps"),
                  " Created at: ",attr(res, "created"),
                  " Version: ",   attr(res, "version"))
    mtext(text=txt, side=1, outer=TRUE, line=-0.5,
          adj=0.02, las=1, cex=0.66)

    par(op)
}

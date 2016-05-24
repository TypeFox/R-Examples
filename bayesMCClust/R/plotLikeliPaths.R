plotLikeliPaths <-
function(outList, from=10, by=1 ) {

    dev.new(width=10, height=10) 

    logDensSeq <- seq(from, outList$Mcmc$M, by)

    noOfRows <- length( intersect( 
                    c("logLike", "logEtaPrior", "logBetaPrior", "logXiPrior", "logEPrior", "logPostDens", "logClassLike", "entropy" ), 
                    names(outList) ) )
    
    par(mfrow=c(noOfRows,1), mai=c(0.1, 0.6, 0.3, 0.1), omi=c(0.3, 0.0, 0.0, 0.1)) # c(bottom, left, top, right)

    if ( is.element("logLike", names(outList)) ) { # 
        plot(logDensSeq, outList$logLike[logDensSeq], type="l", ylab="logLike")
        abline(v=which.max(outList$logLike), col="red", lty=2) 
    }

    if ( is.element("logEtaPrior", names(outList)) )  plot(logDensSeq, outList$logEtaPrior[logDensSeq], type="l", ylab="logEtaPrior")
    if ( is.element("logBetaPrior", names(outList)) ) plot(logDensSeq, outList$logBetaPrior[logDensSeq], type="l", ylab="logBetaPrior")
    if ( is.element("logXiPrior", names(outList)) )   plot(logDensSeq, outList$logXiPrior[logDensSeq], type="l", ylab="logXiPrior")
    if ( is.element("logEPrior", names(outList)) )    plot(logDensSeq, outList$logEPrior[logDensSeq], type="l", ylab="logEPrior")

    if ( is.element("logPostDens", names(outList)) ) {
        plot(logDensSeq, outList$logPostDens[logDensSeq], type="l", ylab="logPostDens")
        abline(v= which.max(outList$logPostDens), col="red", lty=2)
    }

    if ( is.element("logClassLike", names(outList)) ) {
            plot(logDensSeq, outList$logClassLike[logDensSeq], type="l", ylab="logClassLike")
            abline(v= which.max(outList$logClassLike), col="red", lty=2) 
    }

    if ( is.element("entropy", names(outList)) ) {
            plot(logDensSeq, outList$entropy[logDensSeq], type="l", ylab="entropy")
            abline(v= which.min(outList$entropy), col="red", lty=2)
    }

}

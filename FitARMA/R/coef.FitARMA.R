"coef.FitARMA" <-
function (object, ...) 
{
    phiHat <- object$phiHat
    thetaHat <- object$thetaHat
    muHat <- object$muHat
    BETA <- c(phiHat,thetaHat,muHat)
    order <- object$order
    p <- order[1]
    q <- order[3]
    sdB <- sqrt(diag(object$covHat))
    sdmean <- sqrt((object$sigsq)/length(object$res))
    sdfactor <- 1
    sdfactor <- sum(c(1,-thetaHat))/sum(c(1,-phiHat))/sum(c(1,-thetaHat))/sum(c(1,-phiHat))
    sdmean<-sdmean*sdfactor^2
    sdB<-c(sdB,sdmean)
    Z <- BETA/sdB
    rn<-"mu"
    if (q>0)
        rn<-c(paste("theta(", 1:q, ")", sep = ""),rn)
    if (p>0)
        rn<-c(paste("phi(", 1:p, ")", sep = ""),rn)
    cn <- c("MLE", "sd", "Z-ratio")
    ans <- matrix(c(BETA, sdB, Z), ncol = 3)
    dimnames(ans) <- list(rn, cn)
    ans
}


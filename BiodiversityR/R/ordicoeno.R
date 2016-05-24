`ordicoeno` <-
function(x, ordiplot, axis=1, legend=FALSE, cex=0.8, ncol=4, ...) {
#    if (!require(mgcv)) {stop("Requires package mgcv")}
    ordiscore <- scores(ordiplot,display="sites")[,axis]
    original <- cbind(x,ordiscore)
    sorted <- original
    seq <- order(ordiscore)
    sorted[1:nrow(original),] <- original[seq,]
    edfs <- array(NA,dim=c(ncol(x)))
    names(edfs) <- colnames(x)
    grDevices::palette(grDevices::rainbow(ncol(x)))
#
    pchtypes <- c(0:ncol(x))
    names(pchtypes) <- pchtypes
    pchtypes <- pchtypes - trunc(pchtypes/26)*26
#
    gammodel <- mgcv::gam(sorted[,1]~s(ordiscore),data=sorted)
    edfs[1] <- summary(gammodel)$edf
    newdata1 <- data.frame(seq(min(sorted$ordiscore), max(sorted$ordiscore), length = 1000))
    newdata2 <- data.frame(seq(min(sorted$ordiscore), max(sorted$ordiscore), length = 20))
    colnames(newdata1) <- colnames(newdata2) <- "ordiscore"
    gamresult1 <- predict(gammodel, newdata1)
    gamresult2 <- predict(gammodel, newdata2)
    graphics::plot(newdata1$ordiscore, gamresult1, type="l", ylim=c(0,max(x)),
        col=1, pch=0, xlab="site score on ordination axis", ylab="species values", ...)
    graphics::points(newdata2$ordiscore, gamresult2, type="p", col=1, pch=pchtypes[1], cex=cex, ...)    
    for (i in 2:ncol(x)) {
        gammodel <- mgcv::gam(sorted[,i]~s(ordiscore), data=sorted)
        gamresult1 <- predict(gammodel, newdata1)
        gamresult2 <- predict(gammodel, newdata2)
        edfs[i] <- summary(gammodel)$edf
        graphics::points(newdata1$ordiscore, gamresult1, type="l", pch=0, col=i, cex=cex, ...)
        graphics::points(newdata2$ordiscore, gamresult2, type="p", pch=pchtypes[i], col=i, cex=cex, ...)
    }
    colnames <- names(edfs)
    edfs <- as.numeric(edfs)
    names(edfs) <- colnames
    if (legend == T) {
        legend("top", legend=colnames, pch=pchtypes[1:ncol(x)], lty=1, col=c(1:ncol(x)), ncol=ncol)
    }
    grDevices::palette("default")
    cat("edfs from GAM models for each species...\n")
    return(edfs)
}


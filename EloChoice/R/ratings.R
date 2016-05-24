# ratings 15_06_21

ratings <- function(x, show="mean", drawplot=TRUE) {
  # 'show' can be one of the following:
  # original - initial sequence
  # mean - mean across all randomizations
  # range - range across all randomizations
  # all - all values
  # var - variance
  # plot - if true, all the info is plotted

  if(drawplot) {
    temp <- x$ratmat
    temp2 <- apply(temp, 2, range, na.rm=T)
    temp2 <- rbind(temp[1,], temp2)
    temp2 <- rbind(colMeans(temp, na.rm=T), temp2)
    temp2 <- temp2[, rev(order(temp2[1,]))]

    plot(0,0, xlim=c(0, ncol(temp)+1), ylim=range(temp), "n", axes=F, xlab="stimulus", ylab="Elo-rating")
    segments(1:ncol(temp2), temp2[3, ], 1:ncol(temp2),temp2[4, ])
    points(1:ncol(temp2), temp2[1, ], pch=16)
    points(1:ncol(temp2), temp2[2, ], pch=16, col="grey", cex=0.8)
    axis(1, at=1:ncol(temp2), labels=colnames(temp2), lwd=NA)
    axis(2, las=1)
    box()
  }

  if(is.null(show)) show <- "donothing"
  if(show=="original") {
    res <- sort(x$ratmat[1,], decreasing = T)
    return(res)
  }

  if(show == "mean") {
    res <- sort(colMeans(x$ratmat, na.rm=T), decreasing = T)
    return(res)
  }

  if(show == "range") {
    #res <- do.call("rbind", lapply(x$mats, function(X)apply(X, 2, function(x)x[!is.na(x)][length(x)-sum(is.na(x))])))
    res <- apply(x$ratmat, 2, range, na.rm=T)
    res <- res[, rev(order(res[1, ]))]
    return(res)
  }

  if(show == "all") {
    temp <- sort(colMeans(x$ratmat, na.rm=T), decreasing = T)
    res <- x$ratmat[, names(temp)]
    return(res)
  }

  if(show == "var") {
    temp <- sort(colMeans(x$ratmat, na.rm=T), decreasing = T)
    res <- apply(x$ratmat, 2, var, na.rm=T)[ names(temp)]
    return(res)
  }
}







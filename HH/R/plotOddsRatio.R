plotOddsRatio <- function(x,
                          ylab="prob(col1 | row1)",
                          xlab="prob(col1 | row2)",
                          alpha=c(.50, .05),
                          col=trellis.par.get("superpose.line")$col,
                          lwd=trellis.par.get("superpose.line")$lwd,
                          lwd.reference=1,
                          ...) {
  if (missing(xlab) && missing(ylab) && length(dimnames(x))==2) {
    col1 <- dimnames(x)[[2]][1]
    row1 <- dimnames(x)[[1]][1]
    row2 <- dimnames(x)[[1]][2]
    ylab <- paste("prob(",col1,"|",row1,")")
    xlab <- paste("prob(",col1,"|",row2,")")
  }
  col <- rep(col, length=length(alpha) + 1)
  tmp <- OddsRatio(x, alpha[1])
  object <-
    xyplot(tmp$prob1 ~ tmp$prob2 + tmp$ci.prob2[,1] + tmp$ci.prob2[,2],
           type="l", lwd=lwd,
           lty=c(1,2,2), col=col[c(1,2,2)],
           xlab=xlab,
           ylab=ylab,
           aspect=1) +
             layer(panel.abline(a=0, b=1, lty=3, lwd=lwd.reference), data=list(lwd.reference=lwd.reference)) +
               layer(panel.points(y=tmp$p1, x=tmp$p2, pch=13, cex=1, col=col[1]),
                     data=list(tmp=tmp, col=col))
  for (i in seq(along=alpha)[-1]) {
    tmp <- OddsRatio(x, alpha[i])
    for (j in 1:2)
      object <- object +
        layer(panel.lines(y=tmp$prob1, x=tmp$ci.prob2[,j],
                          lty=4, col=col[i+1], lwd=lwd), data=list(tmp=tmp, i=i, j=j, col=col))
  }
  update(object,
         legend=list(right=list(
                       fun="draw.key",
                       args=list(key=list(
                                   space="right",
                                   border=TRUE,
                                   text=list(c("  0%", paste(format(100*(1-alpha),0), "%", sep=""), "x=y", "MLE")),
                                   lines=list(
                                     lty=c(1, 2, rep(4,length(alpha)-1), 3),
                                     lwd=c(rep(lwd[1],length(alpha)+1), 1),
                                     col=c(col[1:(length(alpha)+1)], "black", "transparent")),
                                   points=list(
                                     pch=c(rep(0,length(alpha)+2), 13),
                                     col=c(rep("transparent", length(alpha)+2), col[1])),
                                   padding.text=2),
                         draw=FALSE))))
}

## if (FALSE) {
##   data(hypothermia)
##   plotOddsRatioLattice(hypothermia)
##   plotOddsRatioLattice(hypothermia, alpha=.05)
##   plotOddsRatioLattice(hypothermia, alpha=c(.50, .10, .05))
##   plotOddsRatioLattice(hypothermia, alpha=seq(.50, .05, -.05))
##   plotOddsRatioLattice(hypothermia, col="black")
## }

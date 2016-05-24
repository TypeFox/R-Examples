plot.QmethodRes <- function(x, 
                            xlab='z-scores', ylab='statements',
                            pchlist = NULL, colours = NULL,
                            fnames = NULL, legend = TRUE, 
                            dist = TRUE, pchlist.fill = NULL, 
                            leg.pos="bottomright", ...) {
  dfr <- x$zsc
  lowlim <- floor(min(dfr[[1]]))
  highlim <- ceiling(max(dfr))
  if (is.null(pchlist)) {
    pchlist <- c(1, 2, 0, 5, 6, 16, 17, 15, 18, 21, 24, 23, 22, 3, 4, 7, 8, 9)
    pchlist.fill <- c(16, 17, 15, 23, 25, 16, 17, 15, 18, 21, 24, 23, 22, 3, 4, 7, 8, 9)
  }
  if (dist) pts <- qdc.zsc(x)
  nfactors <- length(dfr)
  sta.order <- order(apply(dfr, 1, sd))
  dfr <- dfr[sta.order, ]
  pts <- pts[sta.order, ]
  if (is.null(colours)) colours <- rainbow(length(dfr))
  if (is.null(fnames) & names(x$zsc)[1] == "zsc_f1") fnames <- paste0("Factor ", 1:nfactors)
  if (is.null(fnames) & names(x$zsc)[1] != "zsc_f1") fnames <- names(x$zsc)
  dotchart(dfr[[1]], lcolor=grey(0.4),
           xlim=c(lowlim, highlim),
           ylab=ylab, xlab=xlab, axis=NULL,
           pch=pchlist[[1]], color=colours[[1]], ...)
  for (i in 2:nfactors){
    points(x=dfr[[i]], 1:length(dfr[[i]]), pch = pchlist[i], type = "p", col=colours[[i]], bg=colours[[i]], ...)
  }
  if (dist) {
    for (i in 1:nfactors){
      points(x=pts[,i], 1:length(pts[,i]), pch = pchlist.fill[i], type = "p", col=colours[[i]], bg=colours[[i]], ...)
    }
  }
  axis(side=2, at=1:nrow(dfr), 
       labels=rownames(dfr), 
       las=1, tick=F, line=-0.5, ...)
  abline(v=seq(from=lowlim, to=highlim, by=0.5), col=grey(0.6), lty=3)
  if (legend) {
    legend(leg.pos, 
           legend=fnames, 
           col=colours[1:nfactors], 
           pch=pchlist[1:nfactors], 
           bty="n")
  }
}
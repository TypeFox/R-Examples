 `groupLikelihoodCurves` <-
function(plotT, plotS, plotC, plotD, main=NULL, cex=0.7) {
 plotT <- update(plotT, main=main, sub="d)", par.settings=list(par.ylab.text=list(cex=cex)))
 plotC <- update(plotC, main=main, sub="c)",par.settings=list(par.zlab.text=list(cex=cex)))
 plotD <- update(plotD, main=main, sub="a)",par.settings=list(par.zlab.text=list(cex=cex)))
 plotS <- update(plotS, main=main, sub="b)",par.settings=list(par.zlab.text=list(cex=cex)))
 plot(plotT, split=c(2,2,2,2), more=TRUE)
 plot(plotC, split=c(2,1,2,2), more=TRUE)
 plot(plotS, split=c(1,2,2,2), more=TRUE)
 plot(plotD, split=c(1,1,2,2), more=FALSE)
 }
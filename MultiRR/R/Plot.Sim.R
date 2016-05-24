Plot.Sim <-
function(x){

IDVint1 <- matrix(unlist(x$IDVCV),x$n.sim,4)[,1]
IDVslope1 <- matrix(unlist(x$IDVCV),x$n.sim,4)[,4]
IDCOVinterceptslope1 <- matrix(unlist(x$IDVCV),x$n.sim,4)[,2]

IDVint <- density(IDVint1)
IDVslope <- density(IDVslope1)
IDCOVinterceptslope <- density(IDCOVinterceptslope1)


SeriesVint1 <- matrix(unlist(x$SeriesVCV),x$n.sim,4)[,1]
SeriesVslope1 <- matrix(unlist(x$SeriesVCV),x$n.sim,4)[,4]
SeriesCOVinterceptslope1 <- matrix(unlist(x$SeriesVCV),x$n.sim,4)[,2]

SeriesVint <- density(SeriesVint1)
SeriesVslope <- density(SeriesVslope1)
SeriesCOVinterceptslope <- density(SeriesCOVinterceptslope1)

RepInt1 <- matrix(unlist(x$RInt),x$n.sim,1, byrow=TRUE)
RepSlope1 <-matrix(unlist(x$RSlope),x$n.sim,1, byrow=TRUE)

RepInt <- density(RepInt1)
RepSlope <- density(RepSlope1)


par(mfrow=c(3,3))
plot(IDVint, main=bquote('V'[ind[0]] ==.(unlist(x$SimVCVInd)[1])))
abline(v=unlist(x$SimVCVInd)[1], lty=2)
segments(IDVint1,rep(quantile(IDVint$y, 0.2), length(IDVint1)), IDVint1,rep(0, length(IDVint1)))

plot(SeriesVint, main=bquote('V'[series[0]] ==.(unlist(x$SimVCVSeries)[1])))
abline(v=unlist(x$SimVCVSeries)[1], lty=2)
segments(SeriesVint1,rep(quantile(SeriesVint$y, 0.2), length(SeriesVint1)), SeriesVint1,rep(0, length(SeriesVint1)))


simRepInt <- (unlist(x$SimVCVInd)[1])/((unlist(x$SimVCVSeries)[1])+(unlist(x$SimVCVInd)[1]))
plot(RepInt, main=bquote('R'[intercept] ==.(simRepInt)))
abline(v=simRepInt, lty=2)
segments(RepInt1,rep(quantile(RepInt$y, 0.2), length(RepInt1)), RepInt1,rep(0, length(RepInt1)))



plot(IDVslope, main=bquote('V'[ind[1]] ==.(unlist(x$SimVCVInd)[4])))
abline(v=unlist(x$SimVCVInd)[4], lty=2)
segments(IDVslope1,rep(quantile(IDVslope$y, 0.2), length(IDVslope1)), IDVslope1,rep(0, length(IDVslope1)))




plot(SeriesVslope, main=bquote('V'[series[1]] ==.(unlist(x$SimVCVSeries)[4])))
abline(v=unlist(x$SimVCVSeries)[4], lty=2)
segments(SeriesVslope1,rep(quantile(SeriesVslope$y, 0.2), length(SeriesVslope1)), SeriesVslope1,rep(0, length(SeriesVslope1)))


simRepSlope <- (unlist(x$SimVCVInd)[4])/((unlist(x$SimVCVSeries)[4])+(unlist(x$SimVCVInd)[4]))
plot(RepSlope, main=bquote('R'[slope] ==.(round(simRepSlope,2))))
abline(v=simRepSlope, lty=2)
segments(RepSlope1,rep(quantile(RepSlope$y, 0.2), length(RepSlope1)), RepSlope1,rep(0, length(RepSlope1)))



plot(IDCOVinterceptslope, main=bquote('Cov'[list(ind[0],ind[1])]==.(unlist(x$SimVCVInd)[2])))
abline(v=unlist(x$SimVCVInd)[2], lty=2)
segments(IDCOVinterceptslope1,rep(quantile(IDCOVinterceptslope$y, 0.2), length(IDCOVinterceptslope1)), IDCOVinterceptslope1,rep(0, length(IDCOVinterceptslope1)))


plot(SeriesCOVinterceptslope, main=bquote('Cov'[list(series[0],series[1])]==.(unlist(x$SimVCVSeries)[2])))
abline(v=unlist(x$SimVCVSeries)[2], lty=2)
segments(SeriesCOVinterceptslope1,rep(quantile(SeriesCOVinterceptslope$y, 0.2), length(SeriesCOVinterceptslope1)), SeriesCOVinterceptslope1,rep(0, length(SeriesCOVinterceptslope1)))

}

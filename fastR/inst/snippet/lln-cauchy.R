x <- rcauchy(10000)
runningMean <- cumsum(x) / 1:length(x)
cauchyPlot <- xyplot(runningMean~1:10000,
    ylab="running mean",xlab="n", type="l");
